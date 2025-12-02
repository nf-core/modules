process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a515d04307ea3e0178af75132105cd36c87d0116c6f9daecf81650b973e870fd/data' :
        'community.wave.seqera.io/library/quast:5.3.0--755a216045b6dbdd' }"

    input:
    tuple val(meta) , path(consensus)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gff)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results
    tuple val(meta), path("${prefix}.tsv")               , emit: tsv
    tuple val(meta), path("${prefix}_transcriptome.tsv") , optional: true , emit: transcriptome
    tuple val(meta), path("${prefix}_misassemblies.tsv") , optional: true , emit: misassemblies
    tuple val(meta), path("${prefix}_unaligned.tsv")     , optional: true , emit: unaligned
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ?  "--features $gff" : ''
    def reference = fasta           ?  "-r $fasta"       : ''
    """
    quast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}

    ln -s ${prefix}/report.tsv ${prefix}.tsv
    [ -f  ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ] && ln -s ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv
    [ -f  ${prefix}/contigs_reports/misassemblies_report.tsv         ] && ln -s ${prefix}/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
    [ -f  ${prefix}/contigs_reports/unaligned_report.tsv             ] && ln -s ${prefix}/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | grep "QUAST" | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ? "--features $gff" : ''
    def reference = fasta           ? "-r $fasta" : ''

    """
    mkdir -p $prefix
    touch $prefix/report.tsv
    touch $prefix/report.html
    touch $prefix/report.pdf
    touch $prefix/quast.log
    touch $prefix/transposed_report.txt
    touch $prefix/transposed_report.tex
    touch $prefix/icarus.html
    touch $prefix/report.tex
    touch $prefix/report.txt

    mkdir -p $prefix/basic_stats
    touch $prefix/basic_stats/cumulative_plot.pdf
    touch $prefix/basic_stats/Nx_plot.pdf
    touch $prefix/basic_stats/genome_GC_content_plot.pdf
    touch $prefix/basic_stats/GC_content_plot.pdf

    mkdir -p $prefix/icarus_viewers
    touch $prefix/icarus_viewers/contig_size_viewer.html

    ln -s $prefix/report.tsv ${prefix}.tsv

    if [ $fasta ]; then
        touch $prefix/basic_stats/NGx_plot.pdf
        touch $prefix/basic_stats/gc.icarus.txt

        mkdir -p $prefix/aligned_stats
        touch $prefix/aligned_stats/NAx_plot.pdf
        touch $prefix/aligned_stats/NGAx_plot.pdf
        touch $prefix/aligned_stats/cumulative_plot.pdf

        mkdir -p $prefix/contigs_reports
        touch $prefix/contigs_reports/all_alignments_transcriptome.tsv
        touch $prefix/contigs_reports/contigs_report_transcriptome.mis_contigs.info
        touch $prefix/contigs_reports/contigs_report_transcriptome.stderr
        touch $prefix/contigs_reports/contigs_report_transcriptome.stdout
        touch $prefix/contigs_reports/contigs_report_transcriptome.unaligned.info
        mkdir -p $prefix/contigs_reports/minimap_output
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords.filtered
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords_tmp
        touch $prefix/contigs_reports/minimap_output/transcriptome.sf
        touch $prefix/contigs_reports/minimap_output/transcriptome.unaligned
        touch $prefix/contigs_reports/minimap_output/transcriptome.used_snps
        touch $prefix/contigs_reports/misassemblies_frcurve_plot.pdf
        touch $prefix/contigs_reports/misassemblies_plot.pdf
        touch $prefix/contigs_reports/misassemblies_report.tex
        touch $prefix/contigs_reports/misassemblies_report.tsv
        touch $prefix/contigs_reports/misassemblies_report.txt
        touch $prefix/contigs_reports/transcriptome.mis_contigs.fa
        touch $prefix/contigs_reports/transposed_report_misassemblies.tex
        touch $prefix/contigs_reports/transposed_report_misassemblies.tsv
        touch $prefix/contigs_reports/transposed_report_misassemblies.txt
        touch $prefix/contigs_reports/unaligned_report.tex
        touch $prefix/contigs_reports/unaligned_report.tsv
        touch $prefix/contigs_reports/unaligned_report.txt

        mkdir -p $prefix/genome_stats
        touch $prefix/genome_stats/genome_info.txt
        touch $prefix/genome_stats/transcriptome_gaps.txt
        touch $prefix/icarus_viewers/alignment_viewer.html

        ln -sf ${prefix}/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
        ln -sf ${prefix}/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv
        ln -sf ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv

    fi

    if ([ $fasta ] && [ $gff ]); then
        touch $prefix/genome_stats/features_cumulative_plot.pdf
        touch $prefix/genome_stats/features_frcurve_plot.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | grep "QUAST" | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
