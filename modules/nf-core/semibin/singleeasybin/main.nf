process SEMIBIN_SINGLEEASYBIN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/semibin:2.2.0--pyhdfd78af_0':
        'biocontainers/semibin:2.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("${prefix}/*.csv")              , emit: csv
    tuple val(meta), path("${prefix}/*.h5")               , emit: model, optional: true
    tuple val(meta), path("${prefix}/*.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}/output_bins/*.fa.gz"), emit: output_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    prefix    = task.ext.prefix ?: "${meta.id}"
    """

    SemiBin2 \\
        $args \\
        single_easy_bin \\
        --input-fasta ${fasta} \\
        --input-bam ${bam} \\
        --output ${prefix} \\
        -t $task.cpus \\
        $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin: \$( SemiBin2 --version )
    END_VERSIONS
    """
    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/{contig_bins,recluster_bins_info}.tsv
    touch ${prefix}/{data,data_split}.csv
    mkdir ${prefix}/output_bins
    touch ${prefix}/output_bins/SemiBin_{0,1,2,3}.fa
    gzip  ${prefix}/output_bins/SemiBin*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin: \$( SemiBin2 --version )
    END_VERSIONS
    """
}
