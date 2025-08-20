process SENTIEON_STARALIGN {
    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'


    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1dfe59ef66d7326b43db9ab1f39ce6220b358a311078c949a208f9c9815d4e/data'
        : 'community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5'}"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    val star_ignore_sjdbgtf
    val seq_platform
    val seq_center

    output:
    tuple val(meta), path('*Log.final.out'),                          emit: log_final
    tuple val(meta), path('*Log.out'),                                emit: log_out
    tuple val(meta), path('*Log.progress.out'),                       emit: log_progress
    tuple val(meta), path('*d.out.bam'),                              emit: bam,                optional: true
    tuple val(meta), path("${prefix}.sortedByCoord.out.bam"),         emit: bam_sorted,         optional: true
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam"), emit: bam_sorted_aligned, optional: true
    tuple val(meta), path('*toTranscriptome.out.bam'),                emit: bam_transcript,     optional: true
    tuple val(meta), path('*Aligned.unsort.out.bam'),                 emit: bam_unsorted,       optional: true
    tuple val(meta), path('*fastq.gz'),                               emit: fastq,              optional: true
    tuple val(meta), path('*.tab'),                                   emit: tab,                optional: true
    tuple val(meta), path('*.SJ.out.tab'),                            emit: spl_junc_tab,       optional: true
    tuple val(meta), path('*.ReadsPerGene.out.tab'),                  emit: read_per_gene_tab,  optional: true
    tuple val(meta), path('*.out.junction'),                          emit: junction,           optional: true
    tuple val(meta), path('*.out.sam'),                               emit: sam,                optional: true
    tuple val(meta), path('*.wig'),                                   emit: wig,                optional: true
    tuple val(meta), path('*.bg'),                                    emit: bedgraph,           optional: true
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each { reads1 << it } : reads.eachWithIndex { v, ix -> (ix & 1 ? reads2 : reads1) << v }
    def ignore_gtf = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile ${gtf}"
    def seq_platform_arg = seq_platform ? "'PL:${seq_platform}'" : ""
    def seq_center_arg = seq_center ? "'CN:${seq_center}'" : ""
    attrRG = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:${prefix}' ${seq_center_arg} 'SM:${prefix}' ${seq_platform_arg}"
    def out_sam_type = args.contains('--outSAMtype') ? '' : '--outSAMtype BAM Unsorted'
    mv_unsorted_bam = args.contains('--outSAMtype BAM Unsorted SortedByCoordinate') ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    sentieon STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        ${out_sam_type} \\
        ${ignore_gtf} \\
        ${attrRG} \\
        ${args}

    ${mv_unsorted_bam}

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(sentieon STAR --version | sed -e "s/STAR_//g")
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(sentieon STAR --version | sed -e "s/STAR_//g")
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
