process SEMIBIN_SINGLEEASYBIN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2aa21f74001110a50915b90b72aca51c1e2c804ce45d686e6f4085efa69f8a5b/data'
        : 'community.wave.seqera.io/library/semibin:2.2.1--3214db8e39e5117b'}"

    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("${prefix}/*.csv")                                , emit: csv
    tuple val(meta), path("${prefix}/*.h5")                                 , emit: model           , optional: true
    tuple val(meta), path("${prefix}/*.tsv")                                , emit: tsv
    tuple val(meta), path("${prefix}/output_bins/*.fa.gz")                  , emit: output_fasta
    tuple val("${task.process}"), val('SemiBin'), eval("SemiBin2 --version"), emit: versions_semibin, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    prefix    = task.ext.prefix ?: "${meta.id}"
    """

    SemiBin2 \\
        ${args} \\
        single_easy_bin \\
        --input-fasta ${fasta} \\
        --input-bam ${bam} \\
        --output ${prefix} \\
        -t ${task.cpus} \\
        ${args2}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/{contig_bins,recluster_bins_info}.tsv
    touch ${prefix}/{data,data_split}.csv
    mkdir ${prefix}/output_bins
    touch ${prefix}/output_bins/SemiBin_{0,1,2,3}.fa
    gzip  ${prefix}/output_bins/SemiBin*
    """
}
