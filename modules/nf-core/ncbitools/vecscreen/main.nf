process NCBITOOLS_VECSCREEN {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/biocontainers/ncbi-tools-bin:6.1.20170106-6-deb_cv2"

    input:
    tuple val(meta), path(fasta_file)
    tuple val(adapters_database_meta), path(adapters_database_directory)

    output:
    tuple val(meta), path("${meta.id}.vecscreen.out")    , emit: vecscreen_output
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The VecScreen module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. 1.0 is arbitrarily used here as the version number
    """
    DB=`find -L ${adapters_database_directory} -maxdepth 1 -name "*.nin" | sed 's/\\.nin\$//'`
    vecscreen -d \$DB ${args} -i ${fasta_file} -o ${prefix}.vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vecscreen: 1.0
    END_VERSIONS
    """

    stub:
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. 1.0 is arbitrarily used here as the version number
    """
    touch ${prefix}.vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vecscreen: 1.0
    END_VERSIONS
    """
}
