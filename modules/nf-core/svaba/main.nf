
process SVABA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::svaba=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svaba:1.1.0--h7d7f7ad_2':
        'quay.io/biocontainers/svaba:1.1.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(bwa_index)
    path dbsnp
    path dbsnp_tbi
    path regions

    output:
    tuple val(meta), path("*.svaba*vcf")                    , emit: sv
    tuple val(meta), path("*.svaba*vcf")                    , emit: indel
    tuple val(meta), path("*.svaba.unfiltered*vcf")         , emit: unfiltered_sv
    tuple val(meta), path("*.svaba.unfiltered*vcf")         , emit: unfiltered_indel
    tuple val(meta), path("*.bps.txt.gz")                   , emit: raw_calls
    tuple val(meta), path("*.discordants.txt.gz")           , emit: discordants, optional: true
    tuple val(meta), path("*.log")                          , emit: log
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def bamlist = normalbam ? "-t ${tumorbam} -n ${normalbam}" : "-t ${tumorbam}"
    def dbsnp   = dbsnp ? "-D ${dbsnp}" : ""
    def regions = regions ? "-k ${regions}" : ""
    def bwa     = bwa_index ? "cp -s ${bwa_index}/* ." : ""

    """
    ${bwa}

    svaba \\
        run \\
        $bamlist \\
        -p $task.cpus \\
        $dbsnp \\
        -a $meta.id \\
        -G $fasta \\
        $regions \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(echo \$(svaba --version 2>&1) | sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/' ))
    END_VERSIONS
    """
}
