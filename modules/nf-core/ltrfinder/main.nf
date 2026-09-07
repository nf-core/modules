process LTRFINDER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_finder_parallel:1.1--hdfd78af_0':
        'quay.io/biocontainers/ltr_finder_parallel:1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.scn") , emit: scn
    tuple val(meta), path("*.gff3"), emit: gff
    tuple val("${task.process}"), val("LTR_FINDER_parallel"), eval("LTR_FINDER_parallel -h 2>&1 | sed -n 's/Version: v//p'"), emit: versions_ltr_finder_parallel, topic: versions
    tuple val("${task.process}"), val("ltr_finder"), eval("ltr_finder -h 2>&1 | sed -n 's/ltr_finder v//p'"), emit: versions_ltr_finder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    LTR_FINDER_parallel \\
        -seq ${fasta} \\
        -threads ${task.cpus} \\
        ${args}

    mv "${fasta}.finder.combine.scn"    "${prefix}.scn"
    mv "${fasta}.finder.combine.gff3"   "${prefix}.gff3"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.scn"
    touch "${prefix}.gff3"
    """
}
