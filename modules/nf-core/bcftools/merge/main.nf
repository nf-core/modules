process BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    def container_image = "/bcftools:1.15.1--h0ea216a_0"
    container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(vcfs), path(tbis)
    path bed
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*.{bcf,vcf}{,.gz}"), emit: merged_variants
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    def regions = bed ? "--regions-file $bed" : ""
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"

    """
    bcftools merge \\
        $regions \\
        --threads $task.cpus \\
        --output ${prefix}.${extension} \\
        $args \\
        *.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
