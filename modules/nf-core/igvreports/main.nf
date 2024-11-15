process IGVREPORTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igv-reports:1.12.0--pyh7cba7a3_0':
        'biocontainers/igv-reports:1.12.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(sites), path(tracks), path(tracks_indicies)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.html") , emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = fasta ? "--fasta ${fasta}" : ""
    // If tracks is not null, create a string of the track paths
    def track_arg = tracks ? "--tracks "+ tracks.collect { it.toString() }.join(' ') : ""
    // if "--tracks" is in the args, then add track_string immediately after it in
    // the args string and set the track_arg to ""
    if (args.contains("--tracks") && track_arg) {
        args = args.replace("--tracks", track_arg)
        track_arg = ""
    }

    """
    create_report $sites \
    $args \
    $fasta \
    $track_arg \
    --output ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igvreports: \$(python -c "import igv_reports; print(igv_reports.__version__)")
    END_VERSIONS
    """
}
