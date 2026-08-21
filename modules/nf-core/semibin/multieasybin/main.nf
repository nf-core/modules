process SEMIBIN_MULTIEASYBIN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4c049230e87478349b2ecca53156faedfb077df85a4f6deaa8d0b22bd8a975c4/data'
:         'community.wave.seqera.io/library/semibin:2.4.1--594344a862aee1b8' }"

    input:
    tuple val(meta), path(fasta), path(bam), path(depths)

    output:
    tuple val(meta), path("${prefix}/bins/*.fa.gz")   , emit: bins
    tuple val(meta), path("${prefix}/samples/*/*.csv"), emit: csv
    tuple val(meta), path("${prefix}/samples/*/*.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/SemiBinRun.log") , emit: log
    tuple val("${task.process}"), val('SemiBin'), eval("SemiBin2 --version"), emit: versions_semibin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    prefix    = task.ext.prefix ?: "${meta.id}"

    def input_depth = bam ? "--input-bam ${bam}" : "--abundance ${depths}"
    """
    SemiBin2 \\
        ${args} \\
        multi_easy_bin \\
        --input-fasta ${fasta} \\
        ${input_depth} \\
        --output ${prefix} \\
        -t ${task.cpus} \\
        ${args2}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/samples/S1/
    mkdir -p ${prefix}/samples/S2/
    mkdir -p ${prefix}/bins

    touch ${prefix}/samples/{S1,S2}/{contig_bins,recluster_bins_info}.tsv
    touch ${prefix}/samples/{S1,S2}/{data,data_cov,data_split}.csv
    touch ${prefix}/SemiBinRun.log

    echo "" | gzip > ${prefix}/bins/S1_SemiBin_0.fa.gz
    echo "" | gzip > ${prefix}/bins/S1_SemiBin_1.fa.gz
    echo "" | gzip > ${prefix}/bins/S2_SemiBin_0.fa.gz
    echo "" | gzip > ${prefix}/bins/S2_SemiBin_1.fa.gz
    """
}
