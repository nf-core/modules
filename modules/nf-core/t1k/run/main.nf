process T1K_RUN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c64676623c43855d3a9dd8a5b02a03b13216e454a8bd1e9a7d9d5695daf5030b/data'
        : 'community.wave.seqera.io/library/t1k:1.0.9--793c68c25d680b97'}"

    input:
    tuple val(meta), path(input), path(fasta), path(cordinates), path(barcodes), path(barcodewhitelist)

    output:
    tuple val(meta), path("*_genotype.tsv")       , emit: genotype_tsv
    tuple val(meta), path("*_allele.tsv")         , emit: allele_tsv
    tuple val(meta), path("*_allele.vcf")         , emit: allele_vcf
    tuple val(meta), path("*_candidate*.fq")      , optional: true, emit: candidate_reads
    tuple val(meta), path("*_aligned.fa")         , optional: true, emit: aligned_reads_single
    tuple val(meta), path("*_aligned_{1,2}.fa")   , optional: true, emit: aligned_reads_paired
    tuple val(meta), path("*_barcode_expr.tsv")   , optional: true, emit: barcode_tsv
    tuple val(meta), path("*_candidate_bc.fa")    , optional: true, emit: barcode_candidate
    tuple val(meta), path("*_aligned_bc.fa")      , optional: true, emit: barcode_aligned
    tuple val("${task.process}"), val('t1k'), eval("run-t1k 2>&1 | head -n 1 | cut -d '-' -f 1 | cut -d v -f 2"), topic: versions, emit: versions_t1k

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_bam = input.toString().endsWith('.bam')
    def input_args = ""
    if (is_bam) {
        // BAM input
        input_args = "-b ${input}"
    }
    else if (meta.single_end) {
        // Single-end FASTQ
        input_args = "-u ${input}"
    }
    else if (!meta.single_end && input.size() == 2) {
        // Paired-end FASTQ
        input_args = "-1 ${input[0]} -2 ${input[1]}"
    }
    else if (!meta.single_end && input.size() == 1) {
        // Interleaved FASTQ
        input_args = "-i ${input}"
    }
    else {
        error "Invalid input files for sample ${meta.id}. Please provide either a BAM file, single-end FASTQ, paired-end FASTQ (2 files), or interleaved FASTQ."
    }
    def cordinates_args = cordinates && is_bam ? "-c ${cordinates}" : ''
    def barcodes_args = barcodes ? "--barcode ${barcodes}" : ''
    def barcodewhitelist_args = barcodewhitelist ? "--barcodeWhiteList ${barcodewhitelist} " : ''
    """
    run-t1k \\
        ${args} \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        -f ${fasta} \\
        ${input_args} \\
        ${cordinates_args} \\
        ${barcodes_args} \\
        ${barcodewhitelist_args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}_genotype.tsv
    touch ${prefix}_allele.tsv
    touch ${prefix}_candidate_1.fq
    touch ${prefix}_candidate_2.fq
    touch ${prefix}_aligned_1.fa
    touch ${prefix}_aligned_2.fa
    touch ${prefix}_allele.vcf
    touch ${prefix}_barcode_expr.tsv
    touch ${prefix}_aligned_bc.fa
    touch ${prefix}_candidate_bc.fa
    """
}
