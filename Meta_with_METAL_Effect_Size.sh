# Meta-analysis with METAL Effect Size
# Charles Mordaunt
# 6/15/18

# MARBLES-EARLI TD vs ASD MetaAnalysis effect size
metal

GENOMICCONTROL ON
MARKER Probe
ALLELE RefAllele OtherAllele
EFFECT Effect
STDERR se
SCHEME STDERR
PROCESS METAL_MARBLES_ASD_eff.txt
PROCESS METAL_EARLI_ASD_eff.txt
OUTFILE MetaAnalysis_ASD_eff .tbl
ANALYZE 
QUIT

# MARBLES-EARLI TD vs NonTD MetaAnalysis effect size
metal

GENOMICCONTROL ON
MARKER Probe
ALLELE RefAllele OtherAllele
EFFECT Effect
STDERR se
SCHEME STDERR
PROCESS METAL_MARBLES_NonTD_eff.txt
PROCESS METAL_EARLI_NonTD_eff.txt
OUTFILE MetaAnalysis_NonTD_eff .tbl
ANALYZE 
QUIT