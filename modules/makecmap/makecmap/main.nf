process MAKECMAP_MAKECMAP {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c2'

     if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the makecmap process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/makecmap:${version}"

    input:
    tuple val(meta), path(fasta)
    val enzyme

    output:
    tuple val(meta), path(cmap), emit: cmap
    tuple val(meta), path(cmapkey) , emit: cmapkey
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    fa2cmap_multi_color.pl -i $fasta -e $enzyme 1
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        makecmap: $version
    END_VERSIONS
    """
}
