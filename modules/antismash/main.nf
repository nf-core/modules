process ANTISMASH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::antismash=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash:6.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash:6.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sequence_input)

    output:
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt")                 , optional: true, emit: clusterblast_file
    tuple val(meta), path("$prefix/css/*.css")                               , emit: css_file
    tuple val(meta), path("$prefix/images")                                  , emit: image_directory
    tuple val(meta), path("$prefix/js/*.js")                                 , emit: javascript
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html")   , optional: true, emit: knownclusterblast_html
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt")            , optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg")                 , optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg")            , optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("$prefix/*.gbk")                                   , emit: gbk_input
    tuple val(meta), path("$prefix/*.json")                                  , emit: json_results
    tuple val(meta), path("$prefix/*.log")                                   , emit: log
    tuple val(meta), path("$prefix/*.zip")                                   , emit: zip
    tuple val(meta), path("$prefix/*region*.gbk")                            , emit: gbk_results
    tuple val(meta), path("${prefix}/clusterblastoutput.txt")                , optional: true, emit: clusterblastoutput
    tuple val(meta), path("$prefix/index.html")                              , emit: html
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt")           , optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("$prefix/regions.js")                              , emit: json_sideloading
    path "versions.yml"                                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    ## We specifically do not include annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes
    antismash \\
        $sequence_input \\
        $args \\
        -c $task.cpus \\
        --output-dir $prefix \\
        --genefinding-tool none \\
        --logfile $prefix/${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
