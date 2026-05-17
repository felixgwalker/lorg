# CRISPR Base Editor Window Visualiser

Visualises the activity window of CRISPR base editors along a guide RNA-target duplex, highlighting which nucleotide positions are within the editing window, which bases can be converted, and which flanking bases may be subject to unwanted bystander edits. Base editor window visualisation is a critical step when designing precision edits for de-extinction projects, where installing the exact ancestral nucleotide while avoiding collateral changes to the proxy genome is paramount.

## Inputs

- Guide RNA sequence (20 nt protospacer, 5'→3')
- Target DNA sequence (must include PAM; protospacer + PAM minimum)
- Base editor type identifier from a supported set: `ABE7.10`, `ABE8e`, `BE3`, `BE4max`, `evoAPOBEC`, `AncBE4max`, or a custom window specification
- Parameters: window start position, window end position (1-indexed from PAM-distal end), optional bystander warning threshold

## Outputs

- A colour-coded guide-target duplex diagram showing editable positions, bystander positions, and out-of-window positions (PNG/SVG)
- A position-by-position editability table in CSV format with columns: position, base, in-window flag, editable flag, bystander risk score
- A bystander edit warning list highlighting non-target bases within the activity window
- An edit outcome prediction table listing possible products and their expected frequencies based on published editor efficiency profiles

## Method

Maps the known activity window of the specified base editor (defined as protospacer positions relative to the PAM) onto the guide-target alignment. Identifies all A bases (for ABEs) or C bases (for CBEs) within the window as primary edit targets. Flags additional A or C bases within the window as potential bystander edits. Estimates relative editing efficiency at each position using published positional preference curves for the selected editor variant. Visualises the duplex with a colour scheme distinguishing target, bystander, and silent positions.

## Dependencies

- `matplotlib` — duplex diagram and position plot rendering
- `pandas` — editability and outcome tables
- `biopython` — sequence complement and alignment utilities
- `numpy` — positional efficiency curve calculations
