process SAVANA_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/savana:1.3.7--pyhdfd78af_0'
        : 'quay.io/biocontainers/savana:1.3.7--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(tumour), path(tumour_index), path(normal), path(normal_index)
    tuple val(meta2),path(ref), path(ref_index)

    output:
    tuple val(meta), path("${prefix}.sv_breakpoints.vcf"), emit: sv_breakpoints_vcf
    tuple val(meta), path("${prefix}.sv_breakpoints.bedpe"), emit: sv_breakpoints_bedpe
    tuple val(meta), path("${prefix}.sv_breakpoints_read_support.tsv"), emit: sv_breakpoints_read_support
    tuple val(meta), path("${prefix}.inserted_sequences.fa"), emit: inserted_sequences
    tuple val("${task.process}"), val("savana"), eval("python -c \"import importlib.metadata as m; print(m.version('savana'))\""), emit: versions_savana, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = (task.ext.args ?: '').trim()
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "."
    """
    savana run \\
        --tumour ${tumour} \\
        --normal ${normal} \\
        --ref ${ref} \\
        --ref_index ${ref_index} \\
        --outdir ${outdir} \\
        --sample ${prefix} \\
        --threads ${task.cpus ?: 1} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "."
    """
    mkdir -p ${outdir}
    touch ${outdir}/${prefix}.sv_breakpoints.vcf
    touch ${outdir}/${prefix}.sv_breakpoints.bedpe
    touch ${outdir}/${prefix}.sv_breakpoints_read_support.tsv
    touch ${outdir}/${prefix}.inserted_sequences.fa
    """
}
