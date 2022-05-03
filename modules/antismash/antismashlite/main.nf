process ANTISMASH_ANTISMASHLITE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::antismash-lite=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.0.1--pyhdfd78af_1' :
        'quay.io/biocontainers/antismash-lite:6.0.1--pyhdfd78af_1' }"

    containerOptions {
        workflow.containerEngine == 'singularity' ?
        "-B $antismash_dir:/usr/local/lib/python3.8/site-packages/antismash" :
        workflow.containerEngine == 'docker' ?
        "-v \$PWD/$antismash_dir:/usr/local/lib/python3.8/site-packages/antismash" :
        ''
        }

    input:
    tuple val(meta), path(sequence_input)
    path(databases)
    path(antismash_dir) // Optional input: AntiSMASH installation folder. It is not needed for using this module with conda, but required for docker/singularity (see meta.yml).

    output:
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt")                 , optional: true, emit: clusterblast_file
    tuple val(meta), path("${prefix}/css/*.css")                             , emit: css_file
    tuple val(meta), path("${prefix}/images")                                , emit: image_directory
    tuple val(meta), path("${prefix}/js/*.js")                               , emit: javascript
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html")   , optional: true, emit: knownclusterblast_html
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt")            , optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg")                 , optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg")            , optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("${prefix}/*.gbk")                                 , emit: gbk_input
    tuple val(meta), path("${prefix}/*.json")                                , emit: json_results
    tuple val(meta), path("${prefix}/*.log")                                 , emit: log
    tuple val(meta), path("${prefix}/*.zip")                                 , emit: zip
    tuple val(meta), path("${prefix}/*region*.gbk")                          , emit: gbk_results
    tuple val(meta), path("${prefix}/clusterblastoutput.txt")                , optional: true, emit: clusterblastoutput
    tuple val(meta), path("${prefix}/index.html")                            , emit: html
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt")           , optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("${prefix}/regions.js")                            , emit: json_sideloading
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"

    """
    ## We specifically do not include annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes
    antismash \\
        $args \\
        -c $task.cpus \\
        --output-dir $prefix \\
        --genefinding-tool none \\
        --logfile $prefix/${prefix}.log \\
        --databases $databases \\
        $sequence_input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash-lite: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${prefix}
    touch ${prefix}/clusterblast/stub_c.stub.txt
    touch ${prefix}/css/stub.css
    touch ${prefix}/images
    touch ${prefix}/js/stub.js
    touch ${prefix}/knownclusterblast/regionstub/ctg.stub.html
    touch ${prefix}/knownclusterblast/stub._c.stub.txt
    touch ${prefix}/svg/clusterblast.stub.svg
    touch ${prefix}/svg/knownclusterblast.stub.svg
    touch ${prefix}/stub.gbk
    touch ${prefix}/stub.json
    touch ${prefix}/stub.log
    touch ${prefix}/stub.zip
    touch ${prefix}/stub.region.stub.gbk
    touch ${prefix}/clusterblastoutput.txt
    touch ${prefix}/index.html
    touch ${prefix}/knownclusterblastoutput.txt
    touch ${prefix}/regions.js

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash-lite: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
