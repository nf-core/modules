process COLABFOLD {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/colabfold:1.5.5--pyh7cba7a3_2'
        : 'quay.io/biocontainers/colabfold:1.5.5--pyh7cba7a3_2'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(weights)

    output:
    tuple val(meta), path("${prefix}/*_rank_*.pdb")        , emit: pdb , optional: true
    tuple val(meta), path("${prefix}/*_scores_rank_*.json"), emit: json, optional: true
    tuple val(meta), path("${prefix}/*.a3m")               , emit: a3m , optional: true
    tuple val(meta), path("${prefix}/*.png")               , emit: png , optional: true
    tuple val(meta), path("${prefix}/log.txt")             , emit: log , optional: true
    tuple val("${task.process}"), val('colabfold'), eval("python -c 'import importlib.metadata as m; print(m.version(\"colabfold\"))'"), topic: versions, emit: versions_colabfold

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def data = weights ? "--data ${weights}" : ''
    """
    # Under Docker the container runs as the host UID, which has no home
    # directory; colabfold_batch needs a writable HOME for its cache.
    export HOME="\$PWD"

    if [ ! -s ${fasta} ]; then
        echo "WARN: ${fasta} is empty, no structure predicted" >&2
    else
        colabfold_batch \\
            ${fasta} \\
            ${prefix} \\
            ${data} \\
            ${args}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
    touch ${prefix}/${prefix}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json
    touch ${prefix}/${prefix}.a3m
    touch ${prefix}/${prefix}_plddt.png
    touch ${prefix}/log.txt
    """
}
