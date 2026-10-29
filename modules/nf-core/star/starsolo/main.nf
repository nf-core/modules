// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process STAR_STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b04fee93569be87eedb18990e8d297a7b75f5e6fe1cb768edbc05a1a842e160/data'
:         'community.wave.seqera.io/library/star:2.7.11b--f22f015f7eb1383f' }"

    input:
    tuple val(meta) , val(solotype), path(reads), path(manifest)
    tuple val(meta2), path(opt_whitelist)
    tuple val(meta3), path(index)
    tuple val(meta4), path(gtf)
    tuple val(meta5), path(splice_juntion)

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
    tuple val(meta), path("${prefix}.Aligned.out.bam")                  , emit: bam_unsorted     , optional: true
    tuple val(meta), path("${prefix}*.tab")                             , emit: tab              , optional: true
    tuple val(meta), path("${prefix}.Aligned.out.sam")                  , emit: sam              , optional: true
    tuple val(meta), path("${prefix}.unmapped_{1,2}.fastq.gz")          , emit: fastq            , optional: true
    tuple val(meta), path("${prefix}.Chimeric.out.junction")            , emit: chimeric_junction, optional: true
    tuple val(meta), path('*.wig')                                      , emit: wig              , optional: true
    tuple val(meta), path('*.bg')                                       , emit: bedgraph         , optional: true
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions, emit: versions_star

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_arg = gtf ? "--sjdbGTFfile $gtf" : ""
    def (forward, reverse) = reads.collate(2).transpose()
    zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat" : ""
    def splice_juntion_arg = splice_juntion ? "--sjdbFileChrStartEnd $splice_juntion" : ""
    def whitelist_args = opt_whitelist ? "--soloCBwhitelist ${[opt_whitelist].flatten().join(' ')}" : "--soloCBwhitelist None"
    if (solotype == "SmartSeq") {
        read_files_in = "--readFilesManifest $manifest"
    } else if (solotype in ["CB_UMI_Simple", "CB_UMI_Complex"]) {
        read_files_in = "--readFilesIn ${reverse.join( "," )} ${forward.join( "," )}"
    } else {
        error("Unknown solotype (${solotype})")
    }

    """
    STAR \\
        --genomeDir $index \\
        $read_files_in $zcat \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --soloType $solotype \\
        $whitelist_args $gtf_arg $splice_juntion_arg \\
        $args

    if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi
    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """

    stub:
    prefix           = task.ext.prefix ?: "${meta.id}"
    def args             = task.ext.args ?: ''
    def args_list = (args =~ /[^\s]+/).collect()
    def feature_idx = args_list.indexOf('--soloFeatures')

    def solo_features = feature_idx >= 0
        ? args_list.subList(feature_idx + 1, args_list.findIndexOf(feature_idx + 1) { arg -> arg.startsWith('--') } < 0
            ? args_list.size()
            : args_list.findIndexOf(feature_idx + 1) { arg -> arg.startsWith('--') })
        : ['Gene']
    def solo_cell_filter = args.contains('--soloCellFilter')
        ? args.split('--soloCellFilter')[1].trim().split(/\s+/)[0]
        : 'CellRanger2.2'

    def touch_solo = solo_features.collect { feature ->
        def touch_raw = """
            mkdir -p ${prefix}.Solo.out/${feature}/raw
            echo "" | gzip > ${prefix}.Solo.out/${feature}/raw/matrix.mtx.gz
            echo "" | gzip > ${prefix}.Solo.out/${feature}/raw/barcodes.tsv.gz
            echo "" | gzip > ${prefix}.Solo.out/${feature}/raw/features.tsv.gz
        """

        def touch_filtered = solo_cell_filter != 'None' ? """
            mkdir -p ${prefix}.Solo.out/${feature}/filtered
            echo "" | gzip > ${prefix}.Solo.out/${feature}/filtered/matrix.mtx.gz
            echo "" | gzip > ${prefix}.Solo.out/${feature}/filtered/barcodes.tsv.gz
            echo "" | gzip > ${prefix}.Solo.out/${feature}/filtered/features.tsv.gz
        """ : ''

        touch_raw + touch_filtered
    }.join('\n')

    def touch_summary = solo_features.contains('Gene') ? """
        touch ${prefix}.Solo.out/Gene/Summary.csv
    """ : ''

    def touch_bam        = args.contains('--outSAMtype') && args.contains('BAM') ? """
        ${args.contains('SortedByCoordinate') ? "touch ${prefix}.Aligned.sortedByCoord.out.bam" : ""}
        ${args.contains('Unsorted') ? "touch ${prefix}.Aligned.out.bam" : ""}
    """ : ""
    def touch_transcriptome = args.contains('--quantMode') && args.contains('TranscriptomeSAM') ? "touch ${prefix}.Aligned.toTranscriptome.out.bam" : ""

    """
    mkdir -p ${prefix}.Solo.out
    $touch_solo
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    $touch_bam
    $touch_transcriptome
    $touch_summary
    touch ${prefix}.Aligned.out.sam
    touch ${prefix}.SJ.out.tab
    """
}
