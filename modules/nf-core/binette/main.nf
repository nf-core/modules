process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/binette:1.2.1--cc07d41be4a5b0b2':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(contig2bin), path(contigs)
    tuple val(meta), path(checkm2_db)

    output:
    tuple val(meta), path("final_bins/*.fa.gz")                      , emit: final_bins
    tuple val(meta), path("${prefix}.final_contig_to_bin.tsv")       , emit: contig2bin
    tuple val(meta), path("input_bins_quality_reports/*.tsv")        , emit: input_bins_quality_reports
    tuple val(meta), path("${prefix}.final_bins_quality_reports.tsv"), emit: final_bins_quality_report
    tuple val("${task.process}"), val('binette'), eval("binette --version"), topic: versions, emit: versions_binette

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    binette \\
        --contig2bin_tables ${contig2bin} \\
        --contigs ${contigs} \\
        --outdir . \\
        ${args}

    find final_bins/ -maxdepth 1 -name "*.fa" -type f | while read file; do
        newname="final_bins/${prefix}.\$(basename "\$file")"
        mv "\$file" "\$newname"
        gzip "\$newname"
    done

    find input_bins_quality_reports/ -maxdepth 1 -name "*.tsv" -type f | while read file; do
        newname="final_bins/${prefix}.\$(basename "\$file")"
        mv "\$file" "\$newname"
    done

    mv final_contig_to_bin.tsv ${prefix}/${prefix}.final_contig_to_bin.tsv
    mv final_bins_quality_reports.tsv ${prefix}/${prefix}.final_bins_quality_reports.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p final_bins
    mkdir -p input_bins_quality_reports

    echo "" | gzip > final_bins/${prefix}.binette_bin1.fa.gz
    echo "" | gzip > final_bins/${prefix}.binette_bin2.fa.gz

    touch ${prefix}/${prefix}.final_contig_to_bin.tsv
    touch ${prefix}/${prefix}.final_bins_quality_reports.tsv
    touch ${prefix}/input_bins_quality_reports/input_bins_1.concoct_bins.tsv
    touch ${prefix}/input_bins_quality_reports/input_bins_1.metabat2_bins.tsv
    """
}
