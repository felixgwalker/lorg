'use client'

import { useEffect, useState } from 'react'
import { FeasibilityReport, GRADE_COLORS, GRADE_TEXT, LAYER_LABELS } from '@/lib/types'
import ScoreGauge from '@/components/ScoreGauge'
import LayerRadarChart from '@/components/LayerRadarChart'
import LayerDetailCard from '@/components/LayerDetailCard'
import EthicsPanel from '@/components/EthicsPanel'

export default function Home() {
  const [report, setReport] = useState<FeasibilityReport | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    fetch('/reports/report.json')
      .then((r) => {
        if (!r.ok) throw new Error(`HTTP ${r.status}`)
        return r.json()
      })
      .then(setReport)
      .catch((e) => setError(e.message))
  }, [])

  if (error) {
    return (
      <main className="min-h-screen flex items-center justify-center p-8">
        <div className="rounded-xl border border-red-800 bg-red-900/20 p-6 max-w-md text-center">
          <p className="text-red-400 font-semibold mb-2">Failed to load report</p>
          <p className="text-slate-400 text-sm">{error}</p>
          <p className="text-slate-500 text-xs mt-3">
            Run the Python scorer first to generate a report:
            <br />
            <code className="text-slate-400">python run_scorer.py data/thylacine_case_study/metrics.json --output src/frontend/public/reports/report.json</code>
          </p>
        </div>
      </main>
    )
  }

  if (!report) {
    return (
      <main className="min-h-screen flex items-center justify-center">
        <p className="text-slate-500 animate-pulse">Loading report&hellip;</p>
      </main>
    )
  }

  const { meta, overall, layers } = report
  const scoreColor = GRADE_COLORS[overall.grade] ?? '#94a3b8'
  const scoreText  = GRADE_TEXT[overall.grade]  ?? 'text-slate-400'

  return (
    <main className="min-h-screen bg-slate-950 text-slate-100 px-4 py-8 max-w-5xl mx-auto">

      {/* ── Header ──────────────────────────────────────────────────────── */}
      <header className="mb-8 border-b border-slate-800 pb-6">
        <p className="text-xs text-slate-500 uppercase tracking-widest mb-2">
          Genomic Resurrection Scorer &nbsp;·&nbsp; v{meta.scorer_version}
        </p>
        <h1 className="text-2xl font-bold text-slate-100">
          {meta.common_name_extinct ?? meta.species_extinct}
        </h1>
        <p className="text-slate-400 text-sm mt-1">
          <em>{meta.species_extinct}</em>
          &ensp;&mdash;&ensp;assessed against proxy&ensp;
          <em>{meta.species_proxy}</em>
          {meta.common_name_proxy ? ` (${meta.common_name_proxy})` : ''}
        </p>
        {meta.ancient_dna_source && (
          <p className="text-slate-600 text-xs mt-2">Source: {meta.ancient_dna_source}</p>
        )}
        <p className="text-slate-600 text-xs mt-1">Analysis date: {meta.analysis_date}</p>
      </header>

      {/* ── Feasibility Index + Radar ────────────────────────────────────── */}
      <section className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-8">
        {/* Score Gauge */}
        <div className="rounded-xl border border-slate-800 bg-slate-900/50 p-6 flex flex-col items-center justify-center">
          <p className="text-xs text-slate-500 uppercase tracking-widest mb-4">
            Overall Feasibility Index
          </p>
          <ScoreGauge
            score={overall.feasibility_index}
            grade={overall.grade}
            gradeLabel={overall.grade_label}
          />
          <div className="mt-4 w-full space-y-2">
            {Object.entries(layers).map(([key, layer]) => (
              <div key={key} className="flex items-center gap-2 text-xs">
                <div className="w-2 h-2 rounded-full shrink-0" style={{ backgroundColor: GRADE_COLORS[layer.grade] }} />
                <span className="text-slate-400 flex-1">{LAYER_LABELS[key]}</span>
                <span className="font-mono font-semibold" style={{ color: GRADE_COLORS[layer.grade] }}>
                  {layer.score.toFixed(1)}
                </span>
                <span className="text-slate-600 text-xs w-12 text-right">
                  ×{(layer.weight * 100).toFixed(0)}%
                </span>
              </div>
            ))}
            <div className="border-t border-slate-800 pt-2 flex items-center gap-2 text-xs">
              <div className="w-2 h-2 rounded-full shrink-0" style={{ backgroundColor: scoreColor }} />
              <span className="text-slate-300 flex-1 font-semibold">Feasibility Index</span>
              <span className={`font-mono font-bold text-sm ${scoreText}`}>
                {overall.feasibility_index.toFixed(1)}
              </span>
            </div>
          </div>
        </div>

        {/* Radar chart */}
        <div className="rounded-xl border border-slate-800 bg-slate-900/50 p-6">
          <p className="text-xs text-slate-500 uppercase tracking-widest mb-4">
            Layer Score Profile
          </p>
          <LayerRadarChart report={report} />
        </div>
      </section>

      {/* ── Layer Cards ─────────────────────────────────────────────────── */}
      <section className="mb-8">
        <h2 className="text-xs text-slate-500 uppercase tracking-widest mb-4">
          Assessment Layers
        </h2>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
          {Object.entries(layers).map(([key, layer]) => (
            <LayerDetailCard
              key={key}
              layerKey={key}
              label={LAYER_LABELS[key] ?? key}
              result={layer as any}
            />
          ))}
        </div>
      </section>

      {/* ── Ethics Panel ────────────────────────────────────────────────── */}
      <section className="mb-8">
        <h2 className="text-xs text-slate-500 uppercase tracking-widest mb-4">
          Ethical &amp; Ecological Detail
        </h2>
        <EthicsPanel result={layers.ethical_ecological as any} />
      </section>

      {/* ── Context / divergence stat callouts ──────────────────────────── */}
      {(layers.divergence as any).context && (
        <section className="mb-8">
          <h2 className="text-xs text-slate-500 uppercase tracking-widest mb-4">
            Divergence Statistics
          </h2>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
            {[
              { label: 'Total SNPs',    value: ((layers.divergence as any).context.total_snps as number).toLocaleString() },
              { label: 'Total Indels',  value: ((layers.divergence as any).context.total_indels as number).toLocaleString() },
              { label: 'Structural Variants', value: ((layers.divergence as any).context.total_svs as number).toLocaleString() },
              { label: 'Divergence Time', value: `${(layers.divergence as any).context.mya_divergence} Mya` },
            ].map(({ label, value }) => (
              <div key={label} className="rounded-lg bg-slate-900 border border-slate-800 p-3">
                <p className="text-xs text-slate-500 mb-1">{label}</p>
                <p className="text-base font-mono font-semibold text-slate-200">{value}</p>
              </div>
            ))}
          </div>
        </section>
      )}

      {(layers.edit_burden as any).context && (
        <section className="mb-8">
          <h2 className="text-xs text-slate-500 uppercase tracking-widest mb-4">
            Edit Burden Statistics
          </h2>
          <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
            {[
              { label: 'Min. Functional Edits', value: ((layers.edit_burden as any).context.total_minimum_edits as number).toLocaleString() },
              { label: 'Coding Region Edits',   value: ((layers.edit_burden as any).context.coding_edits as number).toLocaleString() },
              { label: 'Regulatory Edits',      value: ((layers.edit_burden as any).context.regulatory_edits as number).toLocaleString() },
              { label: 'Structural Variants',   value: ((layers.edit_burden as any).context.functional_svs as number).toLocaleString() },
              { label: 'CRISPR Efficiency',     value: `${(((layers.edit_burden as any).context.crispr_efficiency as number) * 100).toFixed(1)}%` },
              { label: 'Years @ Current Rate',  value: `${(layers.edit_burden as any).context.years_at_current_throughput}` },
            ].map(({ label, value }) => (
              <div key={label} className="rounded-lg bg-slate-900 border border-slate-800 p-3">
                <p className="text-xs text-slate-500 mb-1">{label}</p>
                <p className="text-base font-mono font-semibold text-slate-200">{value}</p>
              </div>
            ))}
          </div>
        </section>
      )}

      {/* Footer */}
      <footer className="border-t border-slate-800 pt-6 text-xs text-slate-600 space-y-1">
        <p>Genomic Resurrection Scorer v{meta.scorer_version} &mdash; reproducible, evidence-based de-extinction feasibility assessment.</p>
        <p>Scores reflect published data and documented weighting logic. Not a recommendation for or against any de-extinction programme.</p>
      </footer>
    </main>
  )
}
