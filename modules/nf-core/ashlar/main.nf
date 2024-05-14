process ASHLAR {
    tag '$meta.id'
    label 'process_single'

    conda "bioconda::ashlar=1.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ashlar:1.17.0--pyh5e36f6f_0' :
        'biocontainers/ashlar:1.17.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(images)
    path(opt_dfp)
    path(opt_ffp)

    output:
    tuple val(meta), path("*.ome.tif"), emit: tif
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args           ?: ''
    def prefix        = task.ext.prefix         ?: "${meta.id}"
    def dfp           = opt_dfp                 ? "--dfp $opt_dfp" : ""
    def ffp           = opt_ffp                 ? "--ffp $opt_ffp" : ""
    def num_files     = images instanceof List  ? images.size()    : 1
    def opt_dfp_size  = opt_dfp instanceof List ? opt_dfp.size()   : 1
    def opt_ffp_size  = opt_ffp instanceof List ? opt_ffp.size()   : 1
    def dfp_validated = (opt_dfp_size == 0 || opt_dfp_size == 1 || opt_dfp_size == num_files) ? true : false
    def ffp_validated = (opt_ffp_size == 0 || opt_ffp_size == 1 || opt_ffp_size == num_files) ? true : false

    if ( !dfp_validated ) { error "Please input only zero, one, or N dfp files, where N is the number of input images" }
    if ( !ffp_validated ) { error "Please input only zero, one, or N ffp files, where N is the number of input images" }

    """

    ashlar \\
        -o ${prefix}.ome.tif \\
        $images \\
        $args \\
        $dfp \\
        $ffp

    sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' ${prefix}.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(ashlar --version 2>&1 | sed 's/^.*ashlar //; s/Using.*\$//' )
    END_VERSIONS
    """
}
