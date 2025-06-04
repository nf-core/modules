process ANTISMASH_ANTISMASH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "nf-core/antismash:8.0.1--pyhdfd78af_0"

    input:
    tuple val(meta), path(sequence_input)
    path databases
    path gff

    output:
    tuple val(meta), path("${prefix}/{css,images,js}")                    , emit: html_accessory_files
    tuple val(meta), path("${prefix}/*.gbk")                              , emit: gbk_input
    tuple val(meta), path("${prefix}/*.json")                             , emit: json_results
    tuple val(meta), path("${prefix}/*.log")                              , emit: log
    tuple val(meta), path("${prefix}/*.zip")                              , emit: zip
    tuple val(meta), path("${prefix}/index.html")                         , emit: html
    tuple val(meta), path("${prefix}/regions.js")                         , emit: json_sideloading
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt")              , emit: clusterblast_file          , optional: true
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html"), emit: knownclusterblast_html     , optional: true
    tuple val(meta), path("${prefix}/knownclusterblast/")                 , emit: knownclusterblast_dir      , optional: true
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt")         , emit: knownclusterblast_txt      , optional: true
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg")              , emit: svg_files_clusterblast     , optional: true
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg")         , emit: svg_files_knownclusterblast, optional: true
    tuple val(meta), path("${prefix}/*region*.gbk")                       , emit: gbk_results                , optional: true
    tuple val(meta), path("${prefix}/clusterblastoutput.txt")             , emit: clusterblastoutput         , optional: true
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt")        , emit: knownclusterblastoutput    , optional: true
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    gff_flag = gff ? "--genefinding-gff3 ${gff}" : ""

    """
    ## We specifically do not include on-the-fly annotations (--genefinding-tool none) as
    ## this should be run as a separate module for versioning purposes

    antismash \\
        ${args} \\
        ${gff_flag} \\
        -c ${task.cpus} \\
        --output-dir ${prefix} \\
        --output-basename ${prefix} \\
        --genefinding-tool none \\
        --logfile ${prefix}/${prefix}.log \\
        --databases ${databases} \\
        ${sequence_input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/css
    mkdir ${prefix}/images
    mkdir ${prefix}/js
    touch ${prefix}/NZ_CP069563.1.region001.gbk
    touch ${prefix}/NZ_CP069563.1.region002.gbk
    touch ${prefix}/css/bacteria.css
    touch ${prefix}/genome.gbk
    touch ${prefix}/genome.json
    touch ${prefix}/genome.zip
    touch ${prefix}/images/about.svg
    touch ${prefix}/index.html
    touch ${prefix}/js/antismash.js
    touch ${prefix}/js/jquery.js
    touch ${prefix}/regions.js
    touch ${prefix}/test.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """
}
