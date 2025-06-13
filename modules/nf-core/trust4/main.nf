process TRUST4 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trust4:1.1.5--h5ca1c30_0':
        'biocontainers/trust4:1.1.5--h5ca1c30_0' }"

    input:
    tuple val(meta), path(bam), path(reads)
    path(fasta)
    path(vdj_reference)
    path(barcode_whitelist)
    val(cell_barcode_read)
    val(umi_read)
    val(read_format)
    

    output:
    tuple val(meta), path("*.tsv")                  , emit: tsv
    tuple val(meta), path("*_airr.tsv")             , emit: airr_files
    tuple val(meta), path("${meta.id}_airr.tsv")    , emit: airr_tsv
    tuple val(meta), path("*_report.tsv")           , emit: report_tsv
    tuple val(meta), path("*.fa")                   , emit: fasta
    tuple val(meta), path("*.out")                  , emit: out
    tuple val(meta), path("*.fq")                   , emit: fq
    tuple val(meta), path("**")                     , emit: outs
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_mode = bam ? "-b ${bam}" : ''
    def single_end_mode = reads && meta.single_end ? "-u ${reads}" : ''
    // reference is optional for fastq input
    def reference = vdj_reference ? "--ref ${vdj_reference}" : ""
    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    def paired_end_mode = reads && (meta.single_end == false) ? "-1 ${forward[0]} -2 ${reverse[0]}" : ''
    // read format is optional
    def readFormat = read_format ? "--readFormat ${read_format}" : ''
    // barcodeWhitelist is optional
    def barcodeWhitelist  = barcode_whitelist ? "--barcodeWhitelist ${barcode_whitelist}" : ""
    // add barcode information if present
    if (cell_barcode_read) {
        if (cell_barcode_read == "R1") {
            barcode = "--barcode ${forward[0]}"
        } else if (cell_barcode_read == "R2") {
            barcode = "--barcode ${reverse[0]}"
        }
    } else {
        barcode = ''
    }
    // add umi information if present
    if (umi_read) {
        if (umi_read == "R1") {
            umi = "--UMI ${forward[0]}"
        } else if (umi_read == "R2") {
            umi = "--UMI ${reverse[0]}"
        }
    } else {
        umi = ''
    }

    """
    run-trust4 \\
        ${bam_mode} \\
        ${single_end_mode} \\
        ${paired_end_mode} \\
        ${barcode} \\
        ${readFormat} \\
        ${umi} \\
        -t $task.cpus \\
        -f ${fasta} \\
        -o ${prefix} \\
        ${reference} \\
        ${barcodeWhitelist} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^/TRUST4 using /' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_airr.tsv
    touch ${prefix}_airr_align.tsv
    touch ${prefix}_report.tsv
    touch ${prefix}_assembled_reads.fa
    touch ${prefix}_annot.fa
    touch ${prefix}_cdr3.out
    touch ${prefix}_raw.out
    touch ${prefix}_final.out
    touch ${prefix}_toassemble.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^/TRUST4 using /' )
    END_VERSIONS
    """
}
