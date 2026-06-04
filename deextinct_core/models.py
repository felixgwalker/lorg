"""Core aDNA-aware data model for the de-extinction toolkit."""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any


@dataclass
class DamageProfile:
    """Post-mortem damage characteristics for one aDNA library."""
    sample_id: str
    ct_rate_5prime: list[float] = field(default_factory=list)
    ga_rate_3prime: list[float] = field(default_factory=list)
    mean_fragment_length: float = 0.0
    authenticity_posterior: float = 0.0
    model: str = "geometric"
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "DamageProfile":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class ReconstructedLocus:
    """One reconstructed locus on the extinct target genome."""
    chrom: str
    start: int
    end: int
    sequence: str
    per_site_posteriors: list[float] = field(default_factory=list)
    source_proxy_chrom: str = ""
    source_proxy_start: int = 0
    source_proxy_end: int = 0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ReconstructedLocus":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class TargetReconstruction:
    """Reconstructed extinct-target sequence at one or more loci."""
    target_species: str
    loci: list[ReconstructedLocus] = field(default_factory=list)
    source_samples: list[str] = field(default_factory=list)
    damage_profile: DamageProfile | None = None
    method: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "TargetReconstruction":
        loci = [ReconstructedLocus.from_dict(l) for l in d.get("loci", [])]
        dp = DamageProfile.from_dict(d["damage_profile"]) if d.get("damage_profile") else None
        return cls(
            target_species=d["target_species"],
            loci=loci,
            source_samples=d.get("source_samples", []),
            damage_profile=dp,
            method=d.get("method", ""),
            metadata=d.get("metadata", {}),
        )


@dataclass
class ProxyGenome:
    """The living proxy species genome being engineered."""
    species: str
    assembly_id: str
    fasta_path: str | Path = ""
    annotation_path: str | Path = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["fasta_path"] = str(self.fasta_path)
        d["annotation_path"] = str(self.annotation_path)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "ProxyGenome":
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})
