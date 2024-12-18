process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_2':
        'biocontainers/star:2.7.11b--h43eeafb_2' }"

    input:
    tuple val(meta), val(solotype), path(reads)
    path(opt_whitelist)
    tuple val(meta2), path(index)

    output:
    tuple val(meta),  path("*.Solo.out")         , emit: counts
    tuple val(meta),  path("*Log.final.out")     , emit: log_final
    tuple val(meta),  path("*Log.out")           , emit: log_out
    tuple val(meta),  path("*Log.progress.out")  , emit: log_progress
    tuple val(meta),  path("*/Gene/Summary.csv") , emit: summary
    path "versions.yml"                          , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat": ""

    // Handle solotype argument logic
    switch(solotype) {
        case "CB_UMI_Simple":
            solotype_args = meta.umi_len ? "--soloUMIlen ${meta.umi_len} " : "";
            solotype_args = solotype_args + (opt_whitelist.name != 'NO_FILE' ? "--soloCBwhitelist ${opt_whitelist} " : "--soloCBwhitelist None ");
            solotype_args = solotype_args + (meta.umi_start ? "--soloUMIstart ${meta.umi_start} " : "");
            solotype_args = solotype_args + (meta.cb_len ? "--soloCBlen ${meta.cb_len} " : "");
            solotype_args = solotype_args + (meta.cb_start ? "--soloCBstart ${meta.cb_start} " : "");
            solotype_args = solotype_args + (meta.barcode_len ? "--soloBarcodeReadLength ${meta.barcode_len} " : "");
            solotype_args = solotype_args + (meta.barcode_mate ? "--soloBarcodeMate ${meta.barcode_mate} " : "");
            break
        case "CB_UMI_Complex":
            solotype_args = meta.cb_position ? "--soloCBposition ${meta.cb_position}" : "";
            solotype_args = solotype_args + (opt_whitelist.name != 'NO_FILE' ? "--soloCBwhitelist ${opt_whitelist} " : "--soloCBwhitelist None ");
            solotype_args = solotype_args + (meta.umi_position ? "--soloUMIposition ${meta.umi_position} " : "");
            solotype_args = solotype_args + (meta.adapter_seq ? "--soloAdapterSequence ${meta.adapter_seq} " : "");
            solotype_args = solotype_args + (meta.max_mismatch_adapter ? "--soloAdapterMismatchesNmax ${meta.max_mismatch_adapter} " : "");
            break
        case "SmartSeq":
            solotype_args = "--soloUMIdedup Exact ";
            solotype_args = solotype_args + (meta.strandedness ? "--soloStrand ${meta.strandedness} " : "");
            solotype_args = solotype_args + "--outSAMattrRGline ID:${prefix} ";
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
    mkdir -p ${prefix}.Solo.out/Gene
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.Solo.out/Gene/Summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
