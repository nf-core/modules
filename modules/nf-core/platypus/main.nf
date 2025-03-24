/*
unfortunately need to output the version manually
because platypus CallVariants does not include --version or -v commend
*/
process PLATYPUS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platypus-variant:0.8.1--py27_1':
        'biocontainers/platypus-variant:0.8.1--py27_1' }"

    input:

    tuple val(meta), path(tumor_file), path(tumor_file_bai), path(control_file),  path(control_file_bai)
    path fasta
    path fai
    path skipregions_file

    output:
    tuple val(meta), path('*.vcf.gz')            , emit: vcf
    tuple val(meta), path('*.vcf.gz.tbi')        , emit: tbi
    tuple val(meta), path('*.log')               , emit: log
    path  "versions.yml"                         , emit: version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamlist = control_file ? "${control_file},${tumor_file}" : "${tumor_file}"
    def skipregions = skipregions_file ? "skipRegionsFile=${skipregions_file}" : ""
    def VERSION = '0.8.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    platypus callVariants \\
        --nCPU=${task.cpus}\\
        --bamFiles=$bamlist \\
        --output=${prefix}.vcf \\
        --refFile=$fasta \\
        --logFileName=${prefix}.log \\
        ${skipregions} \\
        $args

    bgzip  --threads ${task.cpus} -c ${prefix}.vcf > ${prefix}.vcf.gz
    tabix  ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platypus: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.8.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.log
    echo | bgzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        platypus: ${VERSION}
    END_VERSIONS
    """
}
