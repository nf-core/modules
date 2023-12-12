process NCBIGENOMEDOWNLOAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh7cba7a3_0' :
        'biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0' }"

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
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}
