process GRIDSS_GRIDSS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gridss=2.13.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta), path(inputs), path(assembly)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(bwa_index)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf, optional:true
    tuple val(meta), path("*.assembly.bam") , emit: assembly, optional:true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def assembly_bam = assembly ? "--assembly ${assembly}" : ""
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""

    """
    ${bwa}

    gridss \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${assembly} \\
        --jvmheap ${task.memory.toGiga() - 1}g \\
        --otherjvmheap ${task.memory.toGiga() - 1}g \\
        ${inputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def steps = args.contains("-s ") ? args.split('-s ')[-1].split(" ")[0] :
                args.contains("--steps ") ? args.split('--steps ')[-1].split(" ")[0] :
                "all"
    def vcf = steps.contains("call") || steps.contains("all") ? "touch ${prefix}.vcf.gz" : ""
    def assembly_bam = steps.contains("assembly") || steps.contains("all") ? "touch ${prefix}.assembly.bam" : ""
    """
    ${vcf}
    ${assembly_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
