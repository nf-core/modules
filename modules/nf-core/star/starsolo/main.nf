process STAR_STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

    input:
    tuple val(meta), val(solotype), path(reads), path(manifest)
    path(opt_whitelist)
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val(star_ignore_sjdbgtf)

    output:
    tuple val(meta), path('*.Solo.out')                  , emit: counts
    tuple val(meta), path("*.Solo.out/Gene*/raw")        , emit: raw_counts
    tuple val(meta), path("*.Solo.out/Gene*/filtered")   , emit: filtered_counts  , optional: true
    tuple val(meta), path("*.Solo.out/Velocyto/raw")     , emit: raw_velocyto     , optional: true
    tuple val(meta), path("*.Solo.out/Velocyto/filtered"), emit: filtered_velocyto, optional: true
    tuple val(meta), path('*Log.final.out')              , emit: log_final
    tuple val(meta), path('*Log.out')                    , emit: log_out
    tuple val(meta), path('*Log.progress.out')           , emit: log_progress
    tuple val(meta), path('*/Gene/Summary.csv')          , emit: summary
    tuple val(meta), path('*d.out.bam')                  , emit: bam              , optional: true
    tuple val(meta), path('*sortedByCoord.out.bam')      , emit: bam_sorted       , optional: true
    tuple val(meta), path('*toTranscriptome.out.bam')    , emit: bam_transcript   , optional: true
    tuple val(meta), path('*Aligned.unsort.out.bam')     , emit: bam_unsorted     , optional: true
    tuple val(meta), path('*fastq.gz')                   , emit: fastq            , optional: true
    tuple val(meta), path('*.tab')                       , emit: tab              , optional: true
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions, emit: versions_star

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ignore_gtf = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"

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
        // Handle single or multiple whitelist arguments
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
        $whitelist_args $ignore_gtf \\
        $args

    if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.Solo.out/Gene
    mkdir -p ${prefix}.Solo.out/Gene/filtered
    mkdir -p ${prefix}.Solo.out/Gene/raw
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.Solo.out/Gene/Summary.csv
    touch ${prefix}.Solo.out/Gene/raw/raw.mtx
    touch ${prefix}.Solo.out/Gene/filtered/filtered.mtx
    touch ${prefix}.sortedByCoord.out.bam
    """
}
