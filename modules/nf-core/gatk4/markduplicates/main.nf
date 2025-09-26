process GATK4_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3d753d93f57969fe76b8628a8dfcd23ef44bccd08c4ced7089c1f94bf47c89f/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel_htslib_samtools:d3becb6465454c35'}"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*cram"),     emit: cram, optional: true
    tuple val(meta), path("*bam"),      emit: bam,  optional: true
    tuple val(meta), path("*.crai"),    emit: crai, optional: true
    tuple val(meta), path("*.bai"),     emit: bai,  optional: true
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bam"
    
    def output_format = prefix.tokenize('.')[-1] == 'cram' ? "cram" : "bam"
    def output_flag = output_format == 'cram' ? "-Ch" : "-bh"

    def output_index = output_format == 'cram' ? "${prefix}.crai" : "${prefix}.bai"

    def input_list = bam.collect { "--INPUT ${it}" }.join(' ')
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def reference2 = fasta ? "-T ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    if (!fasta && output_format == 'cram') error "Fasta reference is required for CRAM output"

    // Using samtools and not Markduplicates to compress to CRAM speeds up computation:
    // https://medium.com/@acarroll.dna/looking-at-trade-offs-in-compression-levels-for-genomics-tools-eec2834e8b94
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        ${input_list} \\
        --COMPRESSION_LEVEL 0 \\
        --OUTPUT /dev/stdout \\
        --METRICS_FILE ${prefix}.metrics \\
        --TMP_DIR . \\
        ${reference} \\
        ${args} \\
        | samtools view ${args2} ${output_flag} ${reference2} -o ${prefix}##idx##${output_index} --write-index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.bam"
    prefix_no_suffix = task.ext.prefix ? prefix.tokenize('.')[0] : "${meta.id}"
    """
    touch ${prefix_no_suffix}.bam
    touch ${prefix_no_suffix}.cram
    touch ${prefix_no_suffix}.cram.crai
    touch ${prefix_no_suffix}.bai
    touch ${prefix}.metrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
