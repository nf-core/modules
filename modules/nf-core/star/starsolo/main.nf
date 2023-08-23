process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0':
        'biocontainers/star:2.7.10b--h9ee0642_0' }"

    input:
    tuple val(meta), val(solotype), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta),  path('*.Solo.out')         , emit: counts
    tuple val(meta),  path('*Log.final.out')     , emit: log_final
    tuple val(meta),  path('*Log.out')           , emit: log_out
    tuple val(meta),  path('*Log.progress.out')  , emit: log_progress
    tuple val(meta),  path('*/Gene/Summary.csv') , emit: summary
    path "versions.yml"                          , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat": ""

    // solotype helper function
    def getParam(param, flag) {
        return meta[param] ? "${flag} ${meta[param]}" : ""
    }

    // Handle solotype argument logic
    switch(solotype) {
        case "CB_UMI_Simple":
            solotype_args = "${getParam('umi_len', '--soloUMIlen')} ";
            solotype_args = "${solotype_args}${getParam('whitelist', '--soloCBwhitelist')} ";
            solotype_args = "${solotype_args}${getParam('umi_start', '--soloUMIstart')} ";
            solotype_args = "${solotype_args}${getParam('cb_len', '--soloCBlen')} ";
            solotype_args = "${solotype_args}${getParam('cb_start', '--soloCBstart')} ";
            solotype_args = "${solotype_args}${getParam('barcode_len', '--soloBarcodeReadLength')} ";
            solotype_args = "${solotype_args}${getParam('barcode_mate', '--soloBarcodeMate')}";
            break
        case "CB_UMI_Complex":
            solotype_args = "${getParam('cb_position', '--soloCBposition')} ";
            solotype_args = "${solotype_args}${getParam('whitelist', '--soloCBwhitelist')} ";
            solotype_args = "${solotype_args}${getParam('umi_position', '--soloUMIposition')} ";
            solotype_args = "${solotype_args}${getParam('adapter_seq', '--soloAdapterSequence')} ";
            solotype_args = "${solotype_args}${getParam('max_mismatch_adapter', '--soloAdapterMismatchesNmax')}";
            break
        case "SmartSeq":
            solotype_args = "--soloUMIdedup Exact ";
            solotype_args = "${solotype_args}${getParam('strandedness', '--soloStrand')} ";
            solotype_args = "${solotype_args}--outSAMattrRGline ID:${prefix}";
            break
        default:
            log.warn("Unknown output solotype (${solotype})");
            break
    }

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

    if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi

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
