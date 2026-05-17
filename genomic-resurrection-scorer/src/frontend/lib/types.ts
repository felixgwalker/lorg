export interface LayerResult {
  score: number
  grade: string
  weight: number
  components: Record<string, number>
  interpretation: string
  flags: string[]
}

export interface DivergenceResult extends LayerResult {
  context: {
    total_snps: number
    total_indels: number
    total_svs: number
    mya_divergence: number
  }
}

export interface EditBurdenResult extends LayerResult {
  context: {
    total_minimum_edits: number
    coding_edits: number
    regulatory_edits: number
    functional_svs: number
    years_at_current_throughput: number
    crispr_efficiency: number
  }
}

export interface EthicalEcologicalResult extends LayerResult {
  welfare_flags: string[]
  regulatory_flags: string[]
  conservation_flags: string[]
  habitat_details: {
    available: boolean
    quality: string
    area_km2: number | null
  }
}

export interface FeasibilityReport {
  meta: {
    species_extinct: string
    species_proxy: string
    common_name_extinct?: string
    common_name_proxy?: string
    ancient_dna_source?: string
    proxy_genome_source?: string
    analysis_date: string
    scorer_version: string
  }
  overall: {
    feasibility_index: number
    grade: string
    grade_label: string
    layer_weights: Record<string, number>
  }
  layers: {
    ancient_dna_quality: LayerResult
    genomic_completeness: LayerResult
    divergence: DivergenceResult
    edit_burden: EditBurdenResult
    ethical_ecological: EthicalEcologicalResult
  }
}

export const GRADE_COLORS: Record<string, string> = {
  A: '#10b981',
  B: '#14b8a6',
  C: '#eab308',
  D: '#f97316',
  F: '#ef4444',
}

export const GRADE_BG: Record<string, string> = {
  A: 'bg-emerald-900/40 border-emerald-700',
  B: 'bg-teal-900/40 border-teal-700',
  C: 'bg-yellow-900/40 border-yellow-700',
  D: 'bg-orange-900/40 border-orange-700',
  F: 'bg-red-900/40 border-red-700',
}

export const GRADE_TEXT: Record<string, string> = {
  A: 'text-emerald-400',
  B: 'text-teal-400',
  C: 'text-yellow-400',
  D: 'text-orange-400',
  F: 'text-red-400',
}

export const LAYER_LABELS: Record<string, string> = {
  ancient_dna_quality:  'Ancient DNA Quality',
  genomic_completeness: 'Genomic Completeness',
  divergence:           'Divergence',
  edit_burden:          'Edit Burden',
  ethical_ecological:   'Ethical / Ecological',
}

export const LAYER_DESCRIPTIONS: Record<string, string> = {
  ancient_dna_quality:  'Read depth, fragment length, coverage breadth, contamination',
  genomic_completeness: 'Fraction of extinct genome recoverable above confidence thresholds',
  divergence:           'Sequence divergence between extinct species and proxy genome',
  edit_burden:          'Number and complexity of edits to convert proxy toward extinct genome',
  ethical_ecological:   'Habitat, ecological role, welfare, regulatory and conservation conflicts',
}
