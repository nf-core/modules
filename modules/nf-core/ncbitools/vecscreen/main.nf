process NCBITOOLS_VECSCREEN {
    tag "$meta.id"
    label 'process_single'

    container "biocontainers/ncbi-tools-bin:v6.1.20170106-6-deb_cv1"
    //container "${ 'biocontainers/ncbi-tools-bin:v6.1.20170106-6-deb_cv1' }"

    input:
    tuple val(meta), path(fasta_file)
    path adapters_database_directory

    output:
    tuple val(meta), path("${meta.id}.vecscreen.out")    , emit: vecscreen_output
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
    DB=`find -L ${adapters_database_directory} -name "*.nin" | sed 's/\\.nin\$//'`
    echo \$DB > db_path.txt
    vecscreen -d \$DB ${args} -i ${fasta_file} -o ${prefix}.vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. Because of this, the name of the container that contains VecScreen is used here to indicate version
    "${task.process}":
        vecscreen: ncbi-tools-bin_v6.1.20170106-6-deb_cv1.img
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.vecscreen.out

    cat <<-END_VERSIONS > versions.yml
    // WARN: VecScreen doesn't output a version number and doesn't appear to have a Github repository. Because of this, the name of the container that contains VecScreen is used here to indicate version
    "${task.process}":
        vecscreen: ncbi-tools-bin_v6.1.20170106-6-deb_cv1.img
    END_VERSIONS
    """
}