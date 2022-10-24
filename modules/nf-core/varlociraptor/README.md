# Varlociraptor

## Workflow

```mermaid

flowchart TD

FASTA[reference fasta] & BAM[input bam] --> ALIGNMENTPROPERTIES[ESTIMATE ALIGNMENT-PROPERTIES] & PREPROCESS


ALIGNMENTPROPERTIES--"alignment-properties.json"--> CALL
PREPROCESS--"bcf"-->CALL

SCENARIO[scenario yml] --> GENERIC
CALL --> GENERIC & TUMORNORMAL --"bcf"--> FILTER_CALLS --> CONTROLFDR & POSTERIORODDS

````
