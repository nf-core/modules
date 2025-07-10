process FUSIONINSPECTOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/139b94a1f10c6e23a8c27eaed1e5a689db978a513d0ee155e74d35f0970814fe/data' :
        'community.wave.seqera.io/library/fusion-inspector_igv-reports_perl-json-xs_pysam_pruned:c6147971d107ab11'}"

    input:
    tuple val(meta), path(reads), path(fusion_list)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*FusionInspector.fusions.tsv")   , emit: tsv
    tuple val(meta), path("fi_workdir/*.gtf")               , emit: out_gtf, optional:true
    tuple val(meta), path("*FusionInspector.log")           , emit: log
    tuple val(meta), path("*html")                          , emit: html
    tuple val(meta), path("*abridged.tsv")                  , emit: abridged_tsv
    tuple val(meta), path("IGV_inputs")                     , emit: igv_inputs
    tuple val(meta), path("fi_workdir")                     , emit: fi_workdir
    tuple val(meta), path("chckpts_dir")                    , emit: chckpts_dir
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta = meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    FusionInspector \\
        --fusions $fusion_list \\
        --genome_lib ${reference} \\
        $fasta \\
        --CPU ${task.cpus} \\
        -O . \\
        --out_prefix $prefix \\
        --vis $args $args2

    # Touch the output files to make sure they exist
    touch FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.abridged.tsv
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.fusion_inspector_web.html
    mkdir -p IGV_inputs
    mkdir -p fi_workdir
    mkdir -p chckpts_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FusionInspector: \$(FusionInspector --version 2>&1 | grep -i 'version' | sed -e 's/FusionInspector version: //' -e 's/[[:space:]]//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch FusionInspector.log
    touch ${prefix}.FusionInspector.fusions.abridged.tsv
    touch ${prefix}.FusionInspector.fusions.tsv
    touch ${prefix}.fusion_inspector_web.html
    mkdir -p chckpts_dir
    touch chckpts_dir/add_FFPM.ok
    touch chckpts_dir/add_splice_info.ok
    touch chckpts_dir/append_microH_info.ok
    touch chckpts_dir/blast_filter.ok
    touch chckpts_dir/coalesce_junc_n_span.ok
    touch chckpts_dir/cp_consol_bam.ok
    touch chckpts_dir/cp_contigs_file_workdir
    touch chckpts_dir/cp_final.ok
    touch chckpts_dir/cp_gtf_file_workdir.ok
    touch chckpts_dir/cp_tracks_json.ok
    touch chckpts_dir/create_fi_igvjs.ok
    touch chckpts_dir/cytoband.ok
    touch chckpts_dir/EM_adj_counts.ok
    touch chckpts_dir/filter_by_frag_threshs.ok
    touch chckpts_dir/final.abridged.ok
    touch chckpts_dir/fusion_annotator.ok
    touch chckpts_dir/fusion_coding_region_effect.ok
    touch chckpts_dir/fusion_contigs.ok
    touch chckpts_dir/fusion_reports_html.ok
    touch chckpts_dir/get_fusion_JUNCTION_reads_from_bam.ok
    touch chckpts_dir/get_fusion_SPANNING_reads_from_bam.ok
    touch chckpts_dir/index_consol_bam.ok
    touch chckpts_dir/init_EM_adj_counts.ok
    touch chckpts_dir/init_spanning_reads_bam.ok
    touch chckpts_dir/mark_dup_reads.ok
    touch chckpts_dir/mark_dups_reads.index.ok
    touch chckpts_dir/merged_contig_fai.ok
    touch chckpts_dir/merged_contig_gtf_to_bed.ok
    touch chckpts_dir/microH.dat.ok
    touch chckpts_dir/prep_igv_extract_junc_reads.ok
    touch chckpts_dir/prep_igv_junc_reads_bam.ok
    touch chckpts_dir/prep_igv_pfam_bed.ok
    touch chckpts_dir/prep_igv_pfam_gff3.ok
    touch chckpts_dir/prep_igv_seqsim_bed.ok
    touch chckpts_dir/prep_igv_seqsim_gff3.ok
    touch chckpts_dir/prep_spanning_reads.ok
    touch chckpts_dir/run_STAR.ok
    touch chckpts_dir/samtools_idx_junc_reads_bam.ok
    touch chckpts_dir/samtools_index_span_reads_bam.ok
    touch chckpts_dir/span_reads_acc.ok
    touch chckpts_dir/${prefix}.bed.bedsort.ok
    touch chckpts_dir/${prefix}.bed.bgzip.ok
    touch chckpts_dir/${prefix}.bed.tabix.ok
    mkdir -p fi_workdir/_STARgenome
    touch fi_workdir/Log.final.out
    touch fi_workdir/Log.out
    touch fi_workdir/Log.progress.out
    touch fi_workdir/microH.dat
    touch fi_workdir/pipeliner.456.cmds
    touch fi_workdir/SJ.out.tab
    touch fi_workdir/star_align.ok
    touch fi_workdir/_STARgenome/exonGeTrInfo.tab
    touch fi_workdir/_STARgenome/exonInfo.tab
    touch fi_workdir/_STARgenome/geneInfo.tab
    touch fi_workdir/_STARgenome/sjdbInfo.txt
    touch fi_workdir/_STARgenome/sjdbList.fromGTF.out.tab
    touch fi_workdir/_STARgenome/sjdbList.out.tab
    touch fi_workdir/_STARgenome/transcriptInfo.tab
    touch fi_workdir/${prefix}.fa
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh.wSpliceInfo
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh.wSpliceInfo.post_blast_filter
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh.wSpliceInfo.post_blast_filter.info
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh.wSpliceInfo.post_blast_filter.post_promisc_filter
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.EMadj.min_frag_thresh.wSpliceInfo.post_blast_filter.post_promisc_filter.info
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.fusion_junction_read_accs
    touch fi_workdir/${prefix}.fusion_preds.coalesced.summary.fusion_spanning_read_accs
    touch fi_workdir/${prefix}.gtf
    touch fi_workdir/${prefix}.igv.Pfam.gff3
    touch fi_workdir/${prefix}.igv.seqsimilar.gff3
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter.EMadj
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter.EMadj.FFPM
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter.EMadj.FFPM.wMicroH
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter.EMadj.FFPM.wMicroH.annotated
    touch fi_workdir/${prefix}.post_blast_and_promiscuity_filter.EMadj.FFPM.wMicroH.annotated.coding_effect
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.bai
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.failed_reads_during_span_analysis
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.fusion_junc_reads.sam
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.fusion_junction_info
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.fusion_spanning_info
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.fusion_span_reads.sam
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.read_align_counts.idx
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.read_align_counts.idx.ok
    touch fi_workdir/${prefix}.star.cSorted.dupsMarked.bam.spanning_reads_want.idx
    touch fi_workdir/${prefix}.star.sortedByCoord.out.bam
    touch fi_workdir/${prefix}.star.sortedByCoord.out.bam.bai
    touch fi_workdir/${prefix}.star.sortedByCoord.out.bam.bai.ok
    touch fi_workdir/${prefix}.star.sortedByCoord.out.bam.ok
    mkdir -p IGV_inputs
    touch IGV_inputs/cytoBand.txt
    touch IGV_inputs/${prefix}.bed
    gzip -c < /dev/null > IGV_inputs/${prefix}.bed.sorted.bed.gz
    touch IGV_inputs/${prefix}.bed.sorted.bed.gz.tbi
    touch IGV_inputs/${prefix}.consolidated.bam
    touch IGV_inputs/${prefix}.consolidated.bam.bai
    touch IGV_inputs/${prefix}.fa
    touch IGV_inputs/${prefix}.fa.fai
    touch IGV_inputs/${prefix}.fusion_inspector_web.json
    touch IGV_inputs/${prefix}.gtf
    touch IGV_inputs/${prefix}.igv.Pfam.bed
    touch IGV_inputs/${prefix}.igv.seqsimilar.bed
    touch IGV_inputs/${prefix}.junction_reads.bam
    touch IGV_inputs/${prefix}.junction_reads.bam.bai
    touch IGV_inputs/${prefix}.ROI.bed
    touch IGV_inputs/${prefix}.spanning_reads.bam
    touch IGV_inputs/${prefix}.spanning_reads.bam.bai
    touch IGV_inputs/tracks.json
    touch IGV_inputs/TrinityFusion.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FusionInspector: \$(FusionInspector --version 2>&1 | grep -i 'version' | sed -e 's/FusionInspector version: //' -e 's/[[:space:]]//g')
    END_VERSIONS
    """
}
