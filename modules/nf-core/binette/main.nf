process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de7fccc12dc09b996ec3b65df6060b4e3ad284088c5491bc3d97c582e8e7c3f6/data':
        'community.wave.seqera.io/library/binette:1.2.1--cc07d41be4a5b0b2' }"

    input:
    tuple val(meta) , path(contig2bin), path(contigs), path(proteins)
    tuple val(meta2), path(checkm2_db)

    output:
    tuple val(meta), path("final_bins/*.fa.gz")                      , emit: final_bins
    tuple val(meta), path("${prefix}.final_contig_to_bin.tsv")       , emit: contig2bin
    tuple val(meta), path("input_bins_quality_reports/*.tsv")        , emit: input_bins_quality_reports
    tuple val(meta), path("${prefix}.final_bins_quality_reports.tsv"), emit: final_bins_quality_report
    tuple val("${task.process}"), val('binette'), eval("binette --version | sed 's/Binette //'"), topic: versions, emit: versions_binette

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_input = proteins ? "--proteins ${proteins}" : ""
    """
    binette \\
        --contig2bin_tables ${contig2bin} \\
        --contigs ${contigs} \\
        ${proteins_input} \\
        --checkm2_db ${checkm2_db} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}

    find final_bins/ -maxdepth 1 -name "*.fa" -type f -exec gzip {} \\;

    find input_bins_quality_reports/ -maxdepth 1 -name "*.tsv" -type f | while read file; do
        newname="input_bins_quality_reports/${prefix}.\$(basename "\$file")"
        mv "\$file" "\$newname"
    done

    mv final_contig_to_bin.tsv ${prefix}.final_contig_to_bin.tsv
    mv final_bins_quality_reports.tsv ${prefix}.final_bins_quality_reports.tsv
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p final_bins
    mkdir -p input_bins_quality_reports

    echo "" | gzip > final_bins/${prefix}_bin1.fa.gz
    echo "" | gzip > final_bins/${prefix}_bin2.fa.gz

    touch ${prefix}.final_contig_to_bin.tsv
    touch ${prefix}.final_bins_quality_reports.tsv
    touch input_bins_quality_reports/input_bins_1.concoct_bins.tsv
    touch input_bins_quality_reports/input_bins_1.metabat2_bins.tsv
    """
}
