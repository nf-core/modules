process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0':
        'biocontainers/star:2.7.10b--h9ee0642_0' }"

    input:
    tuple val(meta), val(solotype), path(reads)
    path index

    output:
    tuple val(meta), path('*.Solo.out')         , emit: counts
    tuple val(meta), path('*Log.final.out')     , emit: log_final
    tuple val(meta), path('*Log.out')           , emit: log_out
    tuple val(meta), path('*Log.progress.out')  , emit: log_progress
    tuple val(meta), path('*/Gene/Summary.csv') , emit: summary
    path "versions.yml"                         , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def umi_len = meta.umi_len ? "--soloUMIlen ${meta.umi_len}" : ""
    def umi_start = meta.umi_start ? "--soloUMIstart ${meta.umi_start}" : ""
    def cb_len = meta.cb_len ? "--soloCBlen ${meta.cb_len}" : ""
    def zcat = reads.any{it.contains(".gz")} ? "--readFilesCommand zcat": ""
    def whitelist = meta.whitelist ? "--soloCBwhitelist ${meta.whitelist}" : "--soloCBwhitelist None"
    def cb_start = meta.cb_start ? "--soloCBstart ${meta.cb_start}" : ""
    def barcode_len = meta.barcode_len ? "--soloBarcodeReadLength ${meta.barcode_len}" : ""
    def barcode_mate = meta.barcode_mate ? "--soloBarcodeMate ${meta.barcode_mate}" : ""
    def cb_position = meta.cb_position ? "--soloCBposition ${meta.cb_position}" : ""
    def umi_position = meta.umi_position ? "--soloUMIposition ${meta.umi_position}" : ""
    def adapter_seq = meta.adapter_seq ? "--soloAdapterSequence ${meta.adapter_seq}" : ""
    def max_mismatch_adapter = meta.max_mismatch_adapter ? "--soloAdapterMismatchesNmax ${meta.max_mismatch_adapter}" : ""
    def strandedness = meta.strandedness ? "--soloStrand ${meta.strandedness}" : ""
    def solotype_args = (solotype == "CB_UMI_Simple") ?
        "${umi_len} ${whitelist} ${umi_start} ${cb_len} ${cb_start} ${barcode_len} ${barcode_mate}" :
        (solotype == "CB_UMI_Complex") ?
        "${cb_position} ${whitelist} ${umi_position} ${adapter_seq} ${max_mismatch_adapter}" :
        (solotype == "SmartSeq") ?
        "--soloUMIdedup Exact ${strandedness} --outSAMattrRGline ID:${prefix}" :
        ""
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reverse.join( "," )} ${forward.join( "," )} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --soloType $solotype \\
        $zcat \\
        $solotype_args \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}.Solo.out/
    touch ${prefix}.Solo.out/Log.final.out
    touch ${prefix}.Solo.out/Log.out
    touch ${prefix}.Solo.out/Log.progress.out
    touch ${prefix}.Solo.out/Summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
