process GATK_UNIFIEDGENOTYPER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk=3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk:3.5--hdfd78af_11':
        'quay.io/biocontainers/gatk:3.5--hdfd78af_11' }"

    input:
    tuple val(meta), path(input), path(index)
    path(fasta)
    path(fai)
    path(dict)
    path(known_vcf)
    path(intervals)
    path(contamination)
    path(dbsnps)
    path(comp)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def contamination_file = contamination ? "-contaminationFile ${contamination}" : ""
    def dbsnps_file = dbsnps ? "--dbsnp ${dbsnps}" : ""
    def comp_file = comp ? "--comp ${comp}" : ""
    def intervals_file = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK RealignerTargetCreator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk3 \\
        -Xmx${avail_mem}g \\
        -nt ${task.cpus} \\
        -T UnifiedGenotyper \\
        -I ${input} \\
        -R ${fasta} \\
        ${contamination_file} \\
        ${dbsnps_file} \\
        ${comp_file} \\
        ${intervals_file}
        -o ${prefix}.vcf \\
        $args

    gzip -n *.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
