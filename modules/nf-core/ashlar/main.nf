process ASHLAR {
    tag '$meta.id'
    label 'process_single'

    conda "bioconda::ashlar=1.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ashlar:1.17.0--pyh5e36f6f_0' :
        'biocontainers/ashlar:1.17.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(file)
    path(opt_dfp)
    path(opt_ffp)

    output:
    tuple val(meta), path("*.ome.tif"), emit: tif
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def dfp = opt_dfp? "--dfp $opt_dfp": ""
    def ffp = opt_ffp? "--ffp $opt_ffp": ""
    def dfp_validated = !(opt_dfp instanceof List) || ((opt_dfp instanceof List) && (opt_dfp.isEmpty() || opt_dfp.size() == 1)) ? true : false
    def ffp_validated = !(opt_ffp instanceof List) || ((opt_ffp instanceof List) && (opt_ffp.isEmpty() || opt_ffp.size() == 1)) ? true : false

    """
    if [ "$dfp_validated" = "false" ]; then
        echo "Error: please input only zero or one dfp files"
        exit 1
    fi

    if [ "$ffp_validated" = "false" ]; then
        echo "Error: please input only zero or one ffp files"
        exit 1
    fi

    ashlar \\
        $file \\
        $args \\
        $dfp \\
        $ffp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ashlar: \$(ashlar --version 2>&1 | sed 's/^.*ashlar //; s/Using.*\$//' )
    END_VERSIONS
    """
}
