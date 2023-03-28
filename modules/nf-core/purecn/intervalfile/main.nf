
process PURECN_INTERVALFILE {
    label 'process_low'

    // TODO: This needs a proper container
    // cf: https://github.com/bioconda/bioconda-recipes/pull/40076
    // cf: https://github.com/BioContainers/multi-package-containers/pull/2554
    conda "bioconda::bioconductor-purecn=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
        path  target_bed
        path  fasta
        val   genome

    output:
        tuple path("*.txt"), emit: baits_intervals
        // Only produced if --export is used
        tuple path("*.bed"), emit: baits_optimized, optional: true
        path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    Rscript /usr/local/lib/R/library/PureCN/extdata/IntervalFile.R --in-file ${target_bed} \\
        --fasta ${fasta} \\
        --out-file baits_intervals.txt \\
        --genome ${genome} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript /usr/local/lib/R/library/PureCN/extdata/PureCN.R --version)
    END_VERSIONS
    """
}
