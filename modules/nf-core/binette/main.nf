process BINETTE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f49d890620286da8d7663f73091c2caa1b7186ae05da2e885ef039f07d628a96/data'
:         'community.wave.seqera.io/library/binette_gzip:3dad8d26ac1fc14c' }"

    input:
    tuple val(meta) , path(contig2bin), path(bindirs), path(contigs), path(proteins)
    tuple val(meta2), path(checkm2_db)

    output:
    tuple val(meta), path("final_bins/*.fa.gz")              , emit: final_bins, optional: true
    tuple val(meta), path("*.final_contig_to_bin.tsv")       , emit: contig2bin
    tuple val(meta), path("input_bins_quality_reports/*.tsv"), emit: input_bins_quality_reports
    tuple val(meta), path("*.final_bins_quality_reports.tsv"), emit: final_bins_quality_report
    tuple val("${task.process}"), val('binette'), eval("binette --version | sed 's/Binette //'"), topic: versions, emit: versions_binette

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if (contig2bin && bindirs) {
        error("Error: contig2bin and bindirs both provided to Binette but they are mutually exclusive! ")
    }

    def input = contig2bin ? "--contig2bin_tables ${contig2bin}" : "--bin-dirs ${bindirs}"
    def proteins_input = proteins ? "--proteins ${proteins}" : ""
    """
    binette \\
        ${input} \\
        --contigs ${contigs} \\
        ${proteins_input} \\
        --checkm2_db ${checkm2_db} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        --outdir . \\
        ${args}

    if [ -d final_bins ]; then
        find final_bins/ -maxdepth 1 -name "*.fa" -type f -exec gzip {} \\;
    fi

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
