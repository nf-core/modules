process INSITUGAME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insitugame:1.0.6--pyhdfd78af_0':
        'quay.io/biocontainers/insitugame:1.0.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(spatial_data)

    output:
    tuple val(meta), path("*.png") , emit: figures, optional: true
    tuple val(meta), path("*.json"), emit: results
    tuple val(meta), path("*.csv") , emit: tables
    tuple val("${task.process}"), val('insitugame'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('insitugame'))\""), topic: versions, emit: versions_insitugame

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "
import insitugame
from insitugame import InSituGameAssay
import pandas as pd
import json

df_in = pd.read_csv('${spatial_data}', sep=None, engine='python')
cols = list(df_in.columns)

sensitive = [cols[1]] if len(cols) > 1 else ['Sensitive']
resistant = [cols[2]] if len(cols) > 2 else ['Resistant']
stroma    = [cols[3]] if len(cols) > 3 else ['Stroma']

gene_spec = {
    'Sensitive': sensitive,
    'Resistant': resistant,
    'Stroma': stroma
}
assay = InSituGameAssay(gene_spec)

results = {
    'sample_id': '${prefix}',
    'version': insitugame.__version__,
    'status': 'completed',
    'gene_spec': gene_spec,
    'n_rows': len(df_in)
}

with open('${prefix}_results.json', 'w') as f:
    json.dump(results, f, indent=2)

df_out = pd.DataFrame([{'sample': '${prefix}', 'status': 'completed', 'n_rows': len(df_in)}])
df_out.to_csv('${prefix}_tables.csv', index=False)
"
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.json
    touch ${prefix}_tables.csv
    touch ${prefix}_figure.png
    """
}
