process ANTISMASH_ANTISMASH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "nf-core/antismash:8.0.0"

    input:
    tuple val(meta), path(sequence_input)
    path databases
    path gff

    output:
    tuple val(meta), path("${prefix}/clusterblast/*_c*.txt"), optional: true, emit: clusterblast_file
    tuple val(meta), path("${prefix}/{css,images,js}"), emit: html_accessory_files
    tuple val(meta), path("${prefix}/knownclusterblast/region*/ctg*.html"), optional: true, emit: knownclusterblast_html
    tuple val(meta), path("${prefix}/knownclusterblast/"), optional: true, emit: knownclusterblast_dir
    tuple val(meta), path("${prefix}/knownclusterblast/*_c*.txt"), optional: true, emit: knownclusterblast_txt
    tuple val(meta), path("${prefix}/svg/clusterblast*.svg"), optional: true, emit: svg_files_clusterblast
    tuple val(meta), path("${prefix}/svg/knownclusterblast*.svg"), optional: true, emit: svg_files_knownclusterblast
    tuple val(meta), path("${prefix}/*.gbk"), emit: gbk_input
    tuple val(meta), path("${prefix}/*.json"), emit: json_results
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val(meta), path("${prefix}/*.zip"), emit: zip
    tuple val(meta), path("${prefix}/*region*.gbk"), optional: true, emit: gbk_results
    tuple val(meta), path("${prefix}/clusterblastoutput.txt"), optional: true, emit: clusterblastoutput
    tuple val(meta), path("${prefix}/index.html"), emit: html
    tuple val(meta), path("${prefix}/knownclusterblastoutput.txt"), optional: true, emit: knownclusterblastoutput
    tuple val(meta), path("${prefix}/regions.js"), emit: json_sideloading
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
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
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
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
