process PURECN_RUN {
    tag "$meta.id"
    label 'process_medium'

    // TODO: This needs a proper container
    // cf: https://github.com/bioconda/bioconda-recipes/pull/40076
    // cf: https://github.com/BioContainers/multi-package-containers/pull/2554
    conda "bioconda::bioconductor-purecn=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(vcf)

    path coverage
    path normal_db
    path intervals

    //TODO: optional inputs
    //path blacklist
    //path mapping_bias

    val genome

    output:
    //TODO: proper output need to be set up
    //tuple val(meta), path("*.txt"), emit: output_purecn
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //TODO: The following script is a minimal test run, for a production pipeline consider to add the followings:
    // --stats-file ${SAMPLEID}_mutect_stats.txt --mapping-bias-file ${mapping_bias} --snp-blacklist ${blacklist}
    // --force --post-optimize --seed 123 --bootstrapn 500 --normal ${normal_coverage} --fun-segmentation PSCBS
    // --cores ${tasks.cpus}
    """
    Rscript PureCN.R \
        --out purecn/output/${meta.id}/ \\
        --tumor ${coverage} \\
        --sampleid ${meta.id} \\
        --vcf ${vcf} \\
        --normaldb ${normal_db} \\
        --intervals ${intervals} \\
        --genome ${genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript /usr/local/lib/R/library/PureCN/extdata/PureCN.R --version)
    END_VERSIONS
    """
}
