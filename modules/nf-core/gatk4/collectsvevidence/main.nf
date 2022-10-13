process GATK4_COLLECTSVEVIDENCE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    def container_image = "/gatk4:4.2.6.1--hdfd78af_0"
                                                          container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(input), path(input_index), path(allele_count_vcf), path(allele_count_vcf_tbi)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.sr.txt.gz")    , emit: split_read_evidence
    tuple val(meta), path("*.sr.txt.gz.tbi"), emit: split_read_evidence_index
    tuple val(meta), path("*.pe.txt.gz")    , emit: paired_end_evidence
    tuple val(meta), path("*.pe.txt.gz.tbi"), emit: paired_end_evidence_index
    tuple val(meta), path("*.ld.txt.gz")    , emit: allele_counts, optional:true
    tuple val(meta), path("*.ld.txt.gz.tbi"), emit: allele_counts_index, optional:true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def allele_vcf = allele_count_vcf ? "--allele-count-file ${prefix}.ld.txt.gz --allele-count-vcf ${allele_count_vcf}" : ""
    def reference  = fasta ? "--reference ${fasta}" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK COLLECTSVEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CollectSVEvidence \\
        ${args} \\
        --input ${input} \\
        --sr-file ${prefix}.sr.txt.gz \\
        --pe-file ${prefix}.pe.txt.gz \\
        ${allele_vcf} \\
        ${reference} \\
        --tmp-dir . \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
