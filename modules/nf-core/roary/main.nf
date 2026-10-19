process ROARY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6e/6eb36cf3874c118d360e8881a61f427cf71cbe7fd3fee5595bf7c2c52996b0d1/data'
        : 'community.wave.seqera.io/library/roary_perl-algorithm-diff_perl-apache-test_perl-app-cpanminus_pruned:09a4b57815a21cac'}"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}/*"), emit: results
    tuple val(meta), path("${prefix}/*.aln"),  emit: aln, optional: true
    tuple val("${task.process}"), val('roary'), eval('roary --version'), emit: versions_roary, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    roary \\
        ${args} \\
        -p ${task.cpus} \\
        -f ${prefix} \\
        ${gff}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/core_gene_alignment.aln
    touch ${prefix}/gene_presence_absence.csv
    touch ${prefix}/gene_presence_absence.Rtab
    touch ${prefix}/accessory_binary_genes.fa
    touch ${prefix}/pan_genome_reference.fa
    touch ${prefix}/summary_statistics.txt
    """
}
