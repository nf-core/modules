process STAR_STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b04fee93569be87eedb18990e8d297a7b75f5e6fe1cb768edbc05a1a842e160/data'
:         'community.wave.seqera.io/library/star:2.7.11b--f22f015f7eb1383f' }"

    input:
    tuple val(meta), val(solotype), path(reads), path(manifest)
    path(opt_whitelist)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}.Solo.out")                         , emit: counts
    tuple val(meta), path("${prefix}.Solo.out/Gene/raw")                , emit: raw_counts
    tuple val(meta), path("${prefix}.Solo.out/Gene/filtered")           , emit: filtered_counts  , optional: true
    tuple val(meta), path("${prefix}.Solo.out/GeneFull/raw")            , emit: raw_genefull     , optional: true
    tuple val(meta), path("${prefix}.Solo.out/GeneFull/filtered")       , emit: filtered_genefull, optional: true
    tuple val(meta), path("${prefix}.Solo.out/Velocyto/raw")            , emit: raw_velocyto     , optional: true
    tuple val(meta), path("${prefix}.Solo.out/Velocyto/filtered")       , emit: filtered_velocyto, optional: true
    tuple val(meta), path("${prefix}.Log.final.out")                    , emit: log_final
    tuple val(meta), path("${prefix}.Log.out")                          , emit: log_out
    tuple val(meta), path("${prefix}.Log.progress.out")                 , emit: log_progress
    tuple val(meta), path("${prefix}.Solo.out/Gene/Summary.csv")        , emit: summary          , optional: true
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam")    , emit: bam_sorted       , optional: true
    tuple val(meta), path("${prefix}.Aligned.toTranscriptome.out.bam")  , emit: bam_transcript   , optional: true
    tuple val(meta), path("${prefix}.Aligned.unsort.out.bam")           , emit: bam_unsorted     , optional: true
    tuple val(meta), path("${prefix}*.tab")                             , emit: tab              , optional: true
    tuple val(meta), path("${prefix}*.out.sam")                         , emit: sam              , optional: true
    tuple val(meta), path('*fastq.gz')                                  , emit: fastq            , optional: true
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions, emit: versions_star

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_arg = gtf ? "--sjdbGTFfile $gtf" : ""

    // Handle solotype argument logic
    if (solotype == "SmartSeq") {
        zcat = ''
        read_files_in = "--readFilesManifest $manifest"
        whitelist_args = opt_whitelist.name != 'NO_FILE' ? "--soloCBwhitelist ${opt_whitelist} " : "--soloCBwhitelist None "
    } else if (solotype == "CB_UMI_Simple") {
        def (forward, reverse) = reads.collate(2).transpose()
        zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat": ""
        read_files_in = "--readFilesIn ${reverse.join( "," )} ${forward.join( "," )}"
        whitelist_args = opt_whitelist.name != 'NO_FILE' ? "--soloCBwhitelist ${opt_whitelist} " : "--soloCBwhitelist None "
    } else if (solotype == "CB_UMI_Complex") {
        def (forward, reverse) = reads.collate(2).transpose()
        zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat": ""
        read_files_in = "--readFilesIn ${reverse.join( "," )} ${forward.join( "," )}"
        def whitelist_arguments = []
        if (opt_whitelist.name != 'NO_FILE') {
            if (opt_whitelist instanceof Path) {
                whitelist_arguments << opt_whitelist
            } else if (opt_whitelist instanceof List) {
                whitelist_arguments.addAll(opt_whitelist)
            }
        }
        whitelist_args = whitelist_arguments ? "--soloCBwhitelist ${whitelist_arguments.join(' ')}" : "--soloCBwhitelist None "
    } else {
        log.warn("Unknown output solotype (${solotype})")
    }

    """
    STAR \\
        --genomeDir $index \\
        $read_files_in $zcat \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --soloType $solotype \\
        $whitelist_args $gtf_arg \\
        $args

    if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.Solo.out/Gene
    mkdir -p ${prefix}.Solo.out/Gene/filtered
    mkdir -p ${prefix}.Solo.out/Gene/raw
    mkdir -p ${prefix}.Solo.out/GeneFull/filtered
    mkdir -p ${prefix}.Solo.out/GeneFull/raw
    mkdir -p ${prefix}.Solo.out/Velocyto/filtered
    mkdir -p ${prefix}.Solo.out/Velocyto/raw
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.Solo.out/Gene/Summary.csv
    echo "" | gzip > ${prefix}.Solo.out/Gene/raw/matrix.mtx.gz
    echo "" | gzip > ${prefix}.Solo.out/Gene/raw/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Gene/raw/features.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Gene/filtered/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Gene/filtered/features.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Gene/filtered/matrix.mtx.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/raw/matrix.mtx.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/raw/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/raw/features.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/filtered/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/filtered/features.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/GeneFull/filtered/matrix.mtx.gz
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Aligned.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/raw/matrix.mtx.gz
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/raw/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/raw/features.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/filtered/matrix.mtx.gz
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/filtered/barcodes.tsv.gz
    echo "" | gzip > ${prefix}.Solo.out/Velocyto/filtered/features.tsv.gz
    """
}
