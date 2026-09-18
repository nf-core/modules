process STARSUITE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://biodepot/star-suite@sha256:0ced4e4941b8a7347c1583598381e6abe25d798df09759c59ede0e6c6dc4e1d3' :
        'docker.io/biodepot/star-suite@sha256:0ced4e4941b8a7347c1583598381e6abe25d798df09759c59ede0e6c6dc4e1d3' }"

    input:
    tuple val(meta),  path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path("${prefix}")                    , emit: quant_results,        optional: true
    tuple val(meta), path('*quant.sf')                    , emit: quant_transcripts,    optional: true
    tuple val(meta), path('*quant.genes.sf')              , emit: quant_genes,          optional: true
    tuple val(meta), path('*quant.genes.tximport.sf')     , emit: quant_genes_tximport, optional: true
    tuple val(meta), path('*Aligned.out.bam')             , emit: bam,                  optional: true
    tuple val(meta), path('*Aligned.sortedByCoord.out.bam'), emit: bam_sorted,           optional: true
    tuple val(meta), path('*fastq.gz')                    , emit: fastq,                optional: true
    tuple val(meta), path('*Log.final.out')               , emit: log_final
    tuple val(meta), path('*Log.out')                     , emit: log_out
    tuple val(meta), path('*Log.progress.out')            , emit: log_progress
    tuple val("${task.process}"), val('starsuite'), eval('STAR --version'), emit: versions_starsuite, topic: versions
    tuple val("${task.process}"), val('star'), eval('STAR --upstream-version'), emit: versions_star, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "STARSUITE_ALIGN does not support Conda. Please use Docker, Singularity, Apptainer, or Podman instead."
    }
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each { read -> reads1 << read } : reads.eachWithIndex { read, read_index -> (read_index & 1 ? reads2 : reads1) << read }
    def include_gtf  = gtf ? "--sjdbGTFfile $gtf" : ''
    def attr_rg      = args.contains('--outSAMattrRGline') ? '' : "--outSAMattrRGline 'ID:$prefix' 'SM:$prefix'"
    def out_sam_type = args.contains('--outSAMtype') ? '' : '--outSAMtype BAM Unsorted'
    def transcriptome_check = args.tokenize().contains('TranscriptVB') ?
        "test -s $index/transcriptome.fa || { echo 'STAR Suite TranscriptVB quantification requires transcriptome.fa in the genome index.' >&2; exit 1; }" : ''
    """
    $transcriptome_check

    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(',')} ${reads2.join(',')} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $include_gtf \\
        $attr_rg \\
        $args

    if [ -f ${prefix}.quant.sf ]; then
        mkdir ${prefix}
        cp ${prefix}.quant.sf ${prefix}/quant.sf
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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "STARSUITE_ALIGN does not support Conda. Please use Docker, Singularity, Apptainer, or Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/quant.sf
    touch ${prefix}.quant.sf
    touch ${prefix}.quant.genes.sf
    touch ${prefix}.quant.genes.tximport.sf
    touch ${prefix}.Aligned.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    """
}
