process UMITOOLS_GROUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::umi_tools=1.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.4--py38hbff2b2d_1' :
        'biocontainers/umi_tools:1.1.4--py38hbff2b2d_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val create_bam
    val get_group_info

    output:
    tuple val(meta), path("*.log")        , emit: log
    tuple val(meta), path("${prefix}.bam"), optional: true, emit: bam
    tuple val(meta), path("*.tsv")        , optional: true, emit: tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def paired  = meta.single_end ? "" : "--paired"
    output_bam  = create_bam      ? "--output-bam -S ${prefix}.bam" : ""
    group_info  = get_group_info  ? "--group-out ${prefix}.tsv"     : ""

    if (create_bam && "$bam" == "${prefix}.bam") { error "Input and output names are the same, set prefix in module configuration to disambiguate!" }

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    PYTHONHASHSEED=0 umi_tools \\
        group \\
        -I $bam \\
        $output_bam \\
        -L ${prefix}.log \\
        $group_info \\
        $paired \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.log
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
