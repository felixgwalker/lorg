'use client'

import { useState } from 'react'
import { LayerResult, GRADE_COLORS, GRADE_BG, GRADE_TEXT, LAYER_DESCRIPTIONS } from '@/lib/types'

interface Props {
  layerKey: string
  label: string
  result: LayerResult
}

function ComponentBar({ name, value }: { name: string; value: number }) {
  const color =
    value >= 80 ? '#10b981' :
    value >= 65 ? '#14b8a6' :
    value >= 50 ? '#eab308' :
    value >= 35 ? '#f97316' : '#ef4444'

  const label = name.replace(/_/g, ' ')

  return (
    <div className="space-y-1">
      <div className="flex justify-between text-xs text-slate-400">
        <span className="capitalize">{label}</span>
        <span style={{ color }}>{value.toFixed(1)}</span>
      </div>
      <div className="h-1.5 rounded-full bg-slate-800 overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-500"
          style={{ width: `${value}%`, backgroundColor: color }}
        />
      </div>
    </div>
  )
}

export default function LayerDetailCard({ layerKey, label, result }: Props) {
  const [expanded, setExpanded] = useState(false)
  const color = GRADE_COLORS[result.grade] ?? '#94a3b8'
  const bgBorder = GRADE_BG[result.grade] ?? 'bg-slate-800/40 border-slate-700'

  return (
    <div className={`rounded-xl border p-4 ${bgBorder} transition-all`}>
      {/* Header */}
      <div className="flex items-start justify-between gap-3">
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2">
            <span
              className="text-xs font-bold px-2 py-0.5 rounded border"
              style={{ color, borderColor: color, backgroundColor: `${color}20` }}
            >
              {result.grade}
            </span>
            <h3 className="text-sm font-semibold text-slate-200 truncate">{label}</h3>
          </div>
          <p className="text-xs text-slate-500 mt-1">{LAYER_DESCRIPTIONS[layerKey]}</p>
        </div>
        <div className="text-right shrink-0">
          <span className="text-2xl font-bold font-mono" style={{ color }}>
            {result.score.toFixed(1)}
          </span>
          <span className="text-slate-500 text-sm">/100</span>
        </div>
      </div>

      {/* Score bar */}
      <div className="mt-3 h-1.5 rounded-full bg-slate-800 overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-700"
          style={{ width: `${result.score}%`, backgroundColor: color }}
        />
      </div>

      {/* Interpretation */}
      <p className="mt-3 text-xs text-slate-400 leading-relaxed">{result.interpretation}</p>

      {/* Flags */}
      {result.flags.length > 0 && (
        <div className="mt-3 flex flex-wrap gap-1.5">
          {result.flags.map((flag) => (
            <span
              key={flag}
              className="text-xs px-2 py-0.5 rounded bg-red-900/30 border border-red-800 text-red-400"
            >
              {flag.replace(/_/g, ' ').toLowerCase()}
            </span>
          ))}
        </div>
      )}

      {/* Component breakdown toggle */}
      {Object.keys(result.components).length > 0 && (
        <button
          onClick={() => setExpanded(!expanded)}
          className="mt-3 text-xs text-slate-500 hover:text-slate-300 transition-colors flex items-center gap-1"
        >
          <span>{expanded ? '▲' : '▼'}</span>
          <span>Component breakdown</span>
        </button>
      )}

      {expanded && (
        <div className="mt-3 space-y-2 border-t border-slate-800 pt-3">
          {Object.entries(result.components).map(([key, val]) => (
            <ComponentBar key={key} name={key} value={val} />
          ))}
        </div>
      )}
    </div>
  )
}
