process PURECN_NORMALDB {
    tag "$meta.id"
    label 'process_low'

    // TODO: This needs a proper container
    // cf: https://github.com/bioconda/bioconda-recipes/pull/40076
    // cf: https://github.com/BioContainers/multi-package-containers/pull/2554
    conda "bioconda::bioconductor-purecn=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path coverage_files
    val  genome

    output:
    // TODO: Proper output needs to be set up
    //tuple val(meta), path("*.bam"), emit: bam

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //TODO: Set the following commands in $args
    // --normal-panel $normal_panel --assay $assay_name
    """
    Rscript NormalDB.R --out-dir normal_dbs \\
        --coverage-files $coverage_files \\
        --genome $genome \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript /usr/local/lib/R/library/PureCN/extdata/PureCN.R --version)
    END_VERSIONS
    """
}
