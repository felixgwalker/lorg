'use client'

import {
  RadarChart,
  Radar,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  ResponsiveContainer,
  Tooltip,
} from 'recharts'
import { FeasibilityReport, GRADE_COLORS } from '@/lib/types'

interface Props {
  report: FeasibilityReport
}

const AXIS_LABELS: Record<string, string> = {
  ancient_dna_quality:  'aDNA Quality',
  genomic_completeness: 'Completeness',
  divergence:           'Divergence',
  edit_burden:          'Edit Burden',
  ethical_ecological:   'Ethics',
}

export default function LayerRadarChart({ report }: Props) {
  const data = Object.entries(report.layers).map(([key, layer]) => ({
    subject: AXIS_LABELS[key] ?? key,
    score: layer.score,
    fullMark: 100,
  }))

  const overallGrade = report.overall.grade
  const fillColor = GRADE_COLORS[overallGrade] ?? '#6366f1'

  return (
    <div className="w-full h-72">
      <ResponsiveContainer width="100%" height="100%">
        <RadarChart data={data} margin={{ top: 10, right: 30, bottom: 10, left: 30 }}>
          <PolarGrid stroke="#1e293b" />
          <PolarAngleAxis
            dataKey="subject"
            tick={{ fill: '#94a3b8', fontSize: 12, fontFamily: 'ui-monospace, monospace' }}
          />
          <PolarRadiusAxis
            angle={90}
            domain={[0, 100]}
            tick={{ fill: '#475569', fontSize: 10 }}
            tickCount={6}
          />
          <Radar
            name="Layer Score"
            dataKey="score"
            stroke={fillColor}
            fill={fillColor}
            fillOpacity={0.25}
            strokeWidth={2}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: '#1e293b',
              border: '1px solid #334155',
              borderRadius: '6px',
              fontFamily: 'ui-monospace, monospace',
              fontSize: '12px',
            }}
            labelStyle={{ color: '#94a3b8' }}
            itemStyle={{ color: '#f1f5f9' }}
            formatter={(value: number) => [`${value.toFixed(1)} / 100`, 'Score']}
          />
        </RadarChart>
      </ResponsiveContainer>
    </div>
  )
}
