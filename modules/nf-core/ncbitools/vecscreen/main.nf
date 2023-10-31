process VECSCREEN {
    tag "$meta.id"
    label 'process_single'

    container "${ 'biocontainers/ncbi-tools-bin:v6.1.20170106-6-deb_cv1' }"

    input:
    tuple val(meta), path(fasta_file)
    val(adapters_database_file)

    output:
    tuple val(meta), path("${prefix}.vecscreen.out")    , emit: vecscreen_output
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The VecScreen module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vecscreen -d ${adapters_database_file} ${args} -i ${fasta_file} -o ${prefix}.vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. Because of this, the name of the container that contains VecScreen is used here to indicate version
    "${task.process}":
        vecscreen: ncbi-tools-bin_v6.1.20170106-6-deb_cv1.img
    END_VERSIONS
    """

    stub:
    """
    touch vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. Because of this, the name of the container that contains VecScreen is used here to indicate version
    "${task.process}":
        vecscreen: ncbi-tools-bin_v6.1.20170106-6-deb_cv1.img
    END_VERSIONS
    """
}