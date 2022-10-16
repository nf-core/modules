process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::ncbi-amrfinderplus=3.10.30" : null)
    def container_image = "ncbi-amrfinderplus:3.10.30--h6e70893_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    output:
    path "amrfinderdb.tar.gz", emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir amrfinderdb
    amrfinder_update -d amrfinderdb
    tar czvf amrfinderdb.tar.gz -C \$(readlink amrfinderdb/latest) ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
