process BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.{bcf,vcf}{,.gz}"), emit: merged_variants
    tuple val(meta), path("*.{tbi,csi}")      , emit: merged_variants_index, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    def regions = bed ? "--regions-file $bed" : ""
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"

    """
    bcftools merge \\
        $args \\
        $regions \\
        --threads $task.cpus \\
        --output ${prefix}.${extension} \\
        $vcfs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
                args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
                args.contains("--write-index") || args.contains("-W") ? "csi" :
                ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
