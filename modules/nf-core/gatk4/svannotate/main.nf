process GATK4_SVANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bed)
    path(fasta)
    path(fasta_fai)
    path(dict)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: annotated_vcf
    tuple val(meta), path("*.vcf.gz.tbi")   , emit: index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def intervals = bed ? "--intervals ${bed}" : ""
    def reference = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK SVAnnotate] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M" SVAnnotate \\
        --variant ${vcf} \\
        --output ${prefix}.vcf.gz \\
        ${intervals} \\
        ${reference} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
