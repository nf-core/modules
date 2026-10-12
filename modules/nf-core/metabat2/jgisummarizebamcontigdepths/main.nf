process METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metabat2:2.17--hd498684_0'
        : 'quay.io/biocontainers/metabat2:2.17--hd498684_0'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt.gz"), emit: depth
    tuple val("${task.process}"), val('metabat2'), eval('metabat2 --help 2>&1 | sed -n "2s/.*:\\([0-9]*\\.[0-9]*\\).*/\\1/p"'), emit: versions_metabat2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OMP_NUM_THREADS=${task.cpus}

    jgi_summarize_bam_contig_depths \\
        --outputDepth ${prefix}.txt \\
        ${args} \\
        ${bam}

    bgzip --threads ${task.cpus} ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.txt.gz
    """
}
