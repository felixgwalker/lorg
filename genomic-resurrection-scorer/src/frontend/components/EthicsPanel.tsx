'use client'

import { EthicalEcologicalResult } from '@/lib/types'

interface Props {
  result: EthicalEcologicalResult
}

const FLAG_STYLES: Record<string, { bg: string; border: string; text: string; prefix: string }> = {
  welfare:      { bg: 'bg-amber-900/30',  border: 'border-amber-700',  text: 'text-amber-400',  prefix: 'Welfare' },
  regulatory:   { bg: 'bg-blue-900/30',   border: 'border-blue-700',   text: 'text-blue-400',   prefix: 'Regulatory' },
  conservation: { bg: 'bg-orange-900/30', border: 'border-orange-700', text: 'text-orange-400', prefix: 'Conservation' },
}

function FlagChip({ label, style }: { label: string; style: typeof FLAG_STYLES[string] }) {
  return (
    <span className={`inline-flex items-center gap-1.5 text-xs px-2.5 py-1 rounded-full border ${style.bg} ${style.border} ${style.text}`}>
      <span className="font-semibold opacity-70">{style.prefix}</span>
      <span>{label.replace(/_/g, ' ')}</span>
    </span>
  )
}

export default function EthicsPanel({ result }: Props) {
  const { habitat_details, welfare_flags, regulatory_flags, conservation_flags } = result

  const habitatQualityColor: Record<string, string> = {
    excellent: 'text-emerald-400',
    good:      'text-teal-400',
    moderate:  'text-yellow-400',
    poor:      'text-orange-400',
    none:      'text-red-400',
  }

  return (
    <div className="rounded-xl border border-slate-800 bg-slate-900/50 p-5 space-y-4">
      <h3 className="text-sm font-semibold text-slate-300 uppercase tracking-widest">
        Ethical &amp; Ecological Detail
      </h3>

      {/* Habitat */}
      <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
        <div className="rounded-lg bg-slate-800/50 p-3">
          <p className="text-xs text-slate-500 mb-1">Habitat Available</p>
          <p className={`text-sm font-semibold ${habitat_details.available ? 'text-emerald-400' : 'text-red-400'}`}>
            {habitat_details.available ? 'Yes' : 'No'}
          </p>
        </div>
        <div className="rounded-lg bg-slate-800/50 p-3">
          <p className="text-xs text-slate-500 mb-1">Habitat Quality</p>
          <p className={`text-sm font-semibold capitalize ${habitatQualityColor[habitat_details.quality] ?? 'text-slate-400'}`}>
            {habitat_details.quality}
          </p>
        </div>
        {habitat_details.area_km2 != null && (
          <div className="rounded-lg bg-slate-800/50 p-3">
            <p className="text-xs text-slate-500 mb-1">Suitable Area</p>
            <p className="text-sm font-semibold text-slate-200">
              {habitat_details.area_km2.toLocaleString()} km&sup2;
            </p>
          </div>
        )}
      </div>

      {/* Flags */}
      {(welfare_flags.length + regulatory_flags.length + conservation_flags.length) > 0 && (
        <div>
          <p className="text-xs text-slate-500 mb-2 uppercase tracking-wide">Active Flags</p>
          <div className="flex flex-wrap gap-2">
            {welfare_flags.map((f) => (
              <FlagChip key={f} label={f} style={FLAG_STYLES.welfare} />
            ))}
            {regulatory_flags.map((f) => (
              <FlagChip key={f} label={f} style={FLAG_STYLES.regulatory} />
            ))}
            {conservation_flags.map((f) => (
              <FlagChip key={f} label={f} style={FLAG_STYLES.conservation} />
            ))}
          </div>
        </div>
      )}

      {welfare_flags.length === 0 && regulatory_flags.length === 0 && conservation_flags.length === 0 && (
        <p className="text-xs text-slate-500 italic">No welfare, regulatory, or conservation flags raised.</p>
      )}
    </div>
  )
}
