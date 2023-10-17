process SURVIVOR_SIMSV {
    tag "simSV"
    label 'process_single'

    conda "bioconda::survivor=1.0.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--hd03093a_2':
        'biocontainers/survivor:1.0.7--hd03093a_2' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(parameters)
    val(snp_mutation_frequency)
    val(sim_reads)

    output:
    tuple val(meta), path("*.txt")          , emit: parameters, optional:true
    tuple val(meta), path("*.vcf")          , emit: vcf, optional:true
    tuple val(meta), path("*.bed")          , emit: bed, optional:true
    tuple val(meta), path("*.fasta")        , emit: fasta, optional:true
    tuple val(meta), path("*.insertions.fa"), emit: insertions, optional:true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id ?: "simSV"

    def create_parameters = parameters ? "" : "SURVIVOR simSV ${prefix}.txt"
    def params_file = parameters ? parameters : "${prefix}.txt"
    def create_vcf = fasta ? "SURVIVOR simSV ${fasta} ${params_file} ${snp_mutation_frequency} ${sim_reads} ${prefix}" : ""

    """
    ${create_parameters}

    ${create_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id ?: "simSV"

    def create_parameters = parameters ? "" : "touch ${prefix}.txt"
    def create_vcf = fasta ? "touch ${prefix}.vcf" : ""

    """
    ${create_parameters}

    ${create_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
