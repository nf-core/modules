process BAMIXCHECKER {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml" //not in conda
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/abenjak/bamixchecker:1.0.1' :
        'docker.io/abenjak/bamixchecker:1.0.1' }"


    input:
    tuple val(meta), path(bams)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta2), path(dict)
    tuple val(meta3), path(target_bed) //optional via non-empty file
    tuple val(meta4), path(nhSNP)  //optional via non-empty file

    output:
    tuple val(meta), path("BAMixChecker/*.pdf")            , emit: heatmap
    tuple val(meta), path("BAMixChecker/Total_result.txt") , emit: total_result
    tuple val(meta), path("BAMixChecker/*.html")           , emit: html
    tuple val(meta), path("BAMixChecker/*_samples.txt")    , emit: results
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    python /BAMixChecker-1.0.1/BAMixChecker.py \\
        -d ./ \\
        -r $fasta \\
        -p $task.cpus \\
        `if [ \$(wc -w ${target_bed} | cut -f1 -d " ") != 0 ]; then echo -b ${target_bed}; fi` \\
        `if [ \$(wc -w ${nhSNP} | cut -f1 -d " ") != 0 ]; then echo -nhSNP ${nhSNP}; fi` \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamixchecker: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p BAMixChecker
    touch BAMixChecker/BAMixChecker_Heatmap.pdf
    touch BAMixChecker/BAMixChecker_Report.html
    touch BAMixChecker/Matched_samples.txt
    touch BAMixChecker/Total_result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamixchecker: $VERSION
    END_VERSIONS
    """
}
