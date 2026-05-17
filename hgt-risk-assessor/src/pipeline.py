"""Top-level pipeline orchestrator."""

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from src import __version__
from src.aggregator import aggregate, SIGNAL_WEIGHTS
from src.config import DEFAULT_WEIGHT_PROFILE
from src.input_parser import parse_input, resolve_host
from src.models import InputFormat, PipelineResult, SignalResult
from src.signals import gc_content, is_proximity, integron, conjugative, prophage

logger = logging.getLogger(__name__)

_SIGNAL_RUNNERS = [
    ("gc_content",   gc_content.run),
    ("is_proximity", is_proximity.run),
    ("integron",     integron.run),
    ("conjugative",  conjugative.run),
    ("prophage",     prophage.run),
]


def run(
    input_path: Path,
    host_id: str,
    input_format: InputFormat = InputFormat.FASTA,
    ref_path: Optional[Path] = None,
    host_gc: Optional[float] = None,
    data_dir: Path = Path("data"),
    output_dir: Path = Path("results"),
    output_prefix: Optional[str] = None,
    threads: int = 4,
    no_network: bool = False,
    entrez_email: str = "",
    entrez_api_key: str = "",
    blast_bin_dir: Optional[Path] = None,
    skip_signals: Optional[list[str]] = None,
    weight_profile: str = DEFAULT_WEIGHT_PROFILE,
    donor_taxon: Optional[str] = None,
    recipient_taxon: Optional[str] = None,
    flanking_sequence: Optional[str] = None,
) -> PipelineResult:
    skip_signals = skip_signals or []

    query = parse_input(input_path, input_format, ref_path)
    logger.info(
        f"Query: {query.identifier}  length={query.length:,} bp  "
        f"GC={query.gc_content:.1%}"
    )

    host = resolve_host(host_id, host_gc, entrez_email, entrez_api_key)
    logger.info(
        f"Host:  {host.identifier}  GC={host.gc_content:.1%}  "
        f"source={host.source}"
    )

    signal_kwargs = dict(
        data_dir=data_dir,
        threads=threads,
        no_network=no_network,
        blast_bin_dir=blast_bin_dir,
    )

    signal_results: list[SignalResult] = []
    for name, runner in _SIGNAL_RUNNERS:
        if name in skip_signals:
            signal_results.append(SignalResult(
                signal_name=name,
                score=None,
                weight=SIGNAL_WEIGHTS[name],
                evidence={},
                warning="Skipped by --skip-signal",
                skipped=True,
            ))
            logger.info(f"Signal {name}: skipped by user")
            continue

        logger.info(f"Signal {name}: running...")
        result = runner(query=query, host=host, **signal_kwargs)
        signal_results.append(result)

        if result.skipped:
            logger.warning(f"Signal {name}: skipped — {result.warning}")
        else:
            logger.info(f"Signal {name}: score={result.score:.3f}")

    aggregation = aggregate(signal_results)
    logger.info(
        f"Risk index: {aggregation.risk_index:.3f} "
        f"({aggregation.risk_level.value})"
    )
    if aggregation.skipped_signals:
        logger.warning(
            f"Skipped signals: {', '.join(aggregation.skipped_signals)} — "
            "risk index may be underestimated."
        )

    # Three-layer model (runs independently of flat signals; never raises)
    three_layer = None
    try:
        from src.scoring.layers import compute_three_layer
        three_layer = compute_three_layer(
            signal_results=signal_results,
            query=query,
            host=host,
            weight_profile_name=weight_profile,
            donor_taxon=donor_taxon,
            recipient_taxon=recipient_taxon,
            flanking_sequence=flanking_sequence,
        )
        logger.info(
            f"Three-layer index: {three_layer.hgt_risk_index:.3f} "
            f"({three_layer.score_band.value})  "
            f"profile={weight_profile}"
        )
    except Exception as exc:
        logger.warning(f"Three-layer model failed (flat result still valid): {exc}")

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    pipeline_result = PipelineResult(
        query=query,
        host=host,
        aggregation=aggregation,
        run_timestamp=timestamp,
        pipeline_version=__version__,
        three_layer=three_layer,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_prefix or f"{query.identifier}_{timestamp}"

    from src.report import generate_html, generate_json
    html_path = output_dir / f"{prefix}_report.html"
    json_path = output_dir / f"{prefix}_report.json"
    generate_html(pipeline_result, html_path)
    generate_json(pipeline_result, json_path)
    logger.info(f"HTML report: {html_path}")
    logger.info(f"JSON report: {json_path}")

    return pipeline_result
