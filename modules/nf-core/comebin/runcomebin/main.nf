process COMEBIN_RUNCOMEBIN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/comebin:1.0.4--hdfd78af_0':
        'biocontainers/comebin:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam, stageAs: "bam/*")

    output:
    tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    tuple val(meta), path("${prefix}/comebin_res.tsv")         , emit: tsv
    tuple val(meta), path("${prefix}/comebin.log")             , emit: log
    tuple val(meta), path("${prefix}/embeddings.tsv")          , emit: embeddings
    tuple val(meta), path("${prefix}/covembeddings.tsv")       , emit: covembeddings
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    prefix                 = task.ext.prefix ?: "${meta.id}"
    def clean_assembly     = assembly.toString() - ~/\.gz$/
    def decompress_contigs = assembly.toString() == clean_assembly ? "" : "gunzip -q -f $assembly"
    def cleanup            = decompress_contigs ? "rm ${clean_assembly}" : ""
    """
    ${decompress_contigs}

    run_comebin.sh \\
        -t ${task.cpus} \\
        -a ${clean_assembly} \\
        -p bam/ \\
        -o . \\
        $args

    mv comebin_res ${prefix}

    find ${prefix}/comebin_res_bins/*.fa -exec gzip {} \\;

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/comebin_res_bins

    echo "" | gzip > ${prefix}/comebin_res_bins/1.fa.gz
    echo "" | gzip > ${prefix}/comebin_res_bins/2.fa.gz

    touch ${prefix}/comebin_res.tsv
    touch ${prefix}/comebin.log
    touch ${prefix}/embeddings.tsv
    touch ${prefix}/covembeddings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """
}
