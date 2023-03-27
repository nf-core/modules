process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::umi_tools=1.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.4--py38hbff2b2d_1' :
        'quay.io/biocontainers/umi_tools:1.1.4--py38hbff2b2d_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val get_output_stats

    output:
    tuple val(meta), path("*.bam")             , emit: bam
    tuple val(meta), path("*edit_distance.tsv"), optional:true, emit: tsv_edit_distance
    tuple val(meta), path("*per_umi.tsv")      , optional:true, emit: tsv_per_umi
    tuple val(meta), path("*per_position.tsv") , optional:true, emit: tsv_umi_per_position
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    def stats = get_output_stats ? "--output-stats $prefix" : ""

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        $stats \\
        $paired \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
