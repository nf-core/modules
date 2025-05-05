process GRIDSS_GENERATEPONBEDPE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta) , path(vcf), path(bedpe), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bwa_index)

    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    tuple val(meta), path("*.bed")  , emit: bed
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_command = vcf ? "INPUT=${vcf}" : ""
    def bedpe_command = bedpe ? "INPUT_BEDPE=${bedpe}" : ""
    def bed_command   = bed ? "INPUT_BED=${bed}" : ""
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""
    def ref = bwa_index ? "REFERENCE_SEQUENCE=${fasta}" : ""
    """
    ${bwa}
    GeneratePonBedpe \\
        ${vcf_command} \\
        ${bedpe_command} \\
        ${bed_command} \\
        ${ref} \\
        OUTPUT_BEDPE=${prefix}.bedpe \\
        OUTPUT_BED=${prefix}.bed \\
        THREADS=$task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GeneratePonBedpe: \$(echo \$(GeneratePonBedpe --version 2>&1))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedpe
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GeneratePonBedpe: \$(echo \$(GeneratePonBedpe --version 2>&1))
    END_VERSIONS
    """
}
