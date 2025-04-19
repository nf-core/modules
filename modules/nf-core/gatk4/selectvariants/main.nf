process GATK4_SELECTVARIANTS {
    tag "$meta.id"
    label 'process_single'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0':
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path fasta
    path fasta_fai
    path dict
    val intervals

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def avail_mem = 3
    if (task.memory) {
        avail_mem = task.memory.toGiga()
    }

    """
    gatk --java-options "-Xmx${avail_mem}g" \\
        SelectVariants \\
        --variant ${vcf} \\
        --reference ${fasta} \\
        --output ${prefix}.vcf.gz \\
        ${interval_command} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: "\$(gatk --version 2>&1 | sed -n 's/^.*(GATK) v\\([^ ]*\\).*/\\1/p')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gatk4: "\$(gatk --version 2>&1 | sed -n 's/^.*(GATK) v\\([^ ]*\\).*/\\1/p')"
  END_VERSIONS
    """
}
