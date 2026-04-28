process NCBIGENOMEDOWNLOAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh7cba7a3_0' :
        'quay.io/biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0' }"

    input:
    val meta
    path accessions
    path taxids
    val groups

    output:
    tuple val(meta), path("*_genomic.gbff.gz")        , emit: gbk     , optional: true
    tuple val(meta), path("*_genomic.fna.gz")         , emit: fna     , optional: true
    tuple val(meta), path("*_rm.out.gz")              , emit: rm      , optional: true
    tuple val(meta), path("*_feature_table.txt.gz")   , emit: features, optional: true
    tuple val(meta), path("*_genomic.gff.gz")         , emit: gff     , optional: true
    tuple val(meta), path("*_protein.faa.gz")         , emit: faa     , optional: true
    tuple val(meta), path("*_protein.gpff.gz")        , emit: gpff    , optional: true
    tuple val(meta), path("*_wgsmaster.gbff.gz")      , emit: wgs_gbk , optional: true
    tuple val(meta), path("*_cds_from_genomic.fna.gz"), emit: cds     , optional: true
    tuple val(meta), path("*_rna.fna.gz")             , emit: rna     , optional: true
    tuple val(meta), path("*_rna_from_genomic.fna.gz"), emit: rna_fna , optional: true
    tuple val(meta), path("*_assembly_report.txt")    , emit: report  , optional: true
    tuple val(meta), path("*_assembly_stats.txt")     , emit: stats   , optional: true
    tuple val("${task.process}"), val('ncbigenomedownload'), eval('ncbi-genome-download --version'), topic: versions, emit: versions_ncbigenomedownload

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def accessions_opt = accessions ? "-A ${accessions}" : ""
    def taxids_opt     = taxids ? "-t ${taxids}" : ""
    """
    ncbi-genome-download \\
        $args \\
        $accessions_opt \\
        $taxids_opt \\
        --output-folder ./ \\
        --flat-output \\
        --parallel $task.cpus \\
        $groups

    """

    stub:
    """
    touch ${meta.id}_genomic.gbff.gz
    touch ${meta.id}_genomic.fna.gz
    touch ${meta.id}_rm.out.gz
    touch ${meta.id}_feature_table.txt.gz
    touch ${meta.id}_genomic.gff.gz
    touch ${meta.id}_protein.faa.gz
    touch ${meta.id}_protein.gpff.gz
    touch ${meta.id}_wgsmaster.gbff.gz
    touch ${meta.id}_cds_from_genomic.fna.gz
    touch ${meta.id}_rna.fna.gz
    touch ${meta.id}_rna_from_genomic.fna.gz
    touch ${meta.id}_assembly_report.txt
    touch ${meta.id}_assembly_stats.txt
    """
}
