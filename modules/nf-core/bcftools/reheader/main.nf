process BCFTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(header), path(samples)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{csi,tbi}")              , emit: index, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fai_argument      = fai ? "--fai $fai" : ""
    def header_argument   = header ? "--header $header" : ""
    def samples_argument  = samples ? "--samples $samples" : ""

    def args2 = task.ext.args2 ?: '--output-type z'
    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    bcftools \\
        reheader \\
        $fai_argument \\
        $header_argument \\
        $samples_argument \\
        $args \\
        --threads $task.cpus \\
        $vcf \\
        | bcftools view \\
        $args2 \\
        --output ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"
    def index = args2.contains("--write-index=tbi") || args2.contains("-W=tbi") ? "tbi" :
                args2.contains("--write-index=csi") || args2.contains("-W=csi") ? "csi" :
                args2.contains("--write-index") || args2.contains("-W") ? "csi" :
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
