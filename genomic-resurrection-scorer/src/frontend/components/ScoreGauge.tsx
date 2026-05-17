'use client'

import { GRADE_COLORS, GRADE_TEXT } from '@/lib/types'

interface Props {
  score: number
  grade: string
  gradeLabel: string
}

export default function ScoreGauge({ score, grade, gradeLabel }: Props) {
  const r = 72
  const cx = 100
  const cy = 105
  const circumference = 2 * Math.PI * r
  // 270° arc, gap at bottom
  const arcLen = circumference * 0.75
  const scoreLen = (score / 100) * arcLen
  const color = GRADE_COLORS[grade] ?? '#94a3b8'

  return (
    <div className="flex flex-col items-center gap-2">
      <svg width="200" height="190" viewBox="0 0 200 190" aria-label={`Feasibility score: ${score}`}>
        {/* Background track — 270° arc, rotated so gap is at bottom */}
        <circle
          cx={cx} cy={cy} r={r}
          fill="none"
          stroke="#1e293b"
          strokeWidth="14"
          strokeDasharray={`${arcLen} ${circumference}`}
          strokeLinecap="round"
          transform={`rotate(135, ${cx}, ${cy})`}
        />
        {/* Score arc */}
        <circle
          cx={cx} cy={cy} r={r}
          fill="none"
          stroke={color}
          strokeWidth="14"
          strokeDasharray={`${scoreLen} ${circumference}`}
          strokeLinecap="round"
          transform={`rotate(135, ${cx}, ${cy})`}
          style={{ transition: 'stroke-dasharray 0.6s ease' }}
        />
        {/* Score label */}
        <text
          x={cx} y={cy - 6}
          textAnchor="middle"
          fill="white"
          fontSize="36"
          fontWeight="700"
          fontFamily="ui-monospace, monospace"
        >
          {score.toFixed(1)}
        </text>
        <text
          x={cx} y={cy + 20}
          textAnchor="middle"
          fill={color}
          fontSize="13"
          fontWeight="600"
          fontFamily="ui-monospace, monospace"
        >
          Grade {grade}
        </text>
        {/* Scale labels */}
        <text x="18" y={cy + 58} textAnchor="middle" fill="#475569" fontSize="11">0</text>
        <text x="182" y={cy + 58} textAnchor="middle" fill="#475569" fontSize="11">100</text>
      </svg>
      <p className={`text-sm font-semibold tracking-wide ${GRADE_TEXT[grade] ?? 'text-slate-400'}`}>
        {gradeLabel}
      </p>
    </div>
  )
}
