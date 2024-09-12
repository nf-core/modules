process GATK4_VARIANTSTOTABLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.5.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.5.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(arguments_file), path(include_intervals), path(exclude_intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.tsv"),        emit: table
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?:''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arguments_file_arg    = arguments_file ? "--arguments_file ${arguments_file}" : ""
    def include_intervals_arg = include_intervals ? "-L ${include_intervals}" : ""
    def exclude_intervals_arg = exclude_intervals ? "-XL ${exclude_intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantsToTable] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantsToTable \\
        ${args} \\
        --variant $vcf \\
        --output ${prefix}.tsv \\
        --reference $fasta \\
        --tmp-dir . \\
        $arguments_file_arg \\
        $include_intervals \\
        $exclude_intervals

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
