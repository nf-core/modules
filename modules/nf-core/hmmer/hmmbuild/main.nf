process HMMER_HMMBUILD {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    def container_image = "hmmer:3.3.2--h87f3376_2"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(alignment)
    path mxfile

    output:
    tuple val(meta), path("*.hmm.gz"), emit: hmm
    path "*.hmmbuild.txt",             emit: hmmbuildout
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def mxfileopt = mxfile ? "--mxfile ${mxfile}" : ""

    """
    hmmbuild \\
        $args \\
        --cpu $task.cpus \\
        -n ${prefix}  \\
        -o ${prefix}.hmmbuild.txt \\
        ${mxfileopt} \\
        ${prefix}.hmm \\
        $alignment

    gzip ${prefix}.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmbuild -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """
}
