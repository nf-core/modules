
process BCFTOOLS_CSQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_1':
        'biocontainers/bcftools:1.20--h8b25389_1' }"


    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(gff3)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.tbi")         , emit: tbi, optional: true
    tuple val(meta), path("*.csi")         , emit: csi, optional: true
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    extension = getVcfExtension(args);

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    bcftools csq \\
        --output ${prefix}.${extension} \\
        --threads ${task.cpus} \\
        --fasta-ref ${fasta} \\
        --gff-annot ${gff3} \\
        $args \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = getVcfExtension(args);
    index = getVcfIndex(args, extension);

    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = index ? "touch ${prefix}.${extension}.${index}" : ""

    if ("$vcf" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
// Custom Functions
String getVcfExtension(String args) {
    return args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
        args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
        args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
        args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
        "vcf";
}
String getVcfIndex(String args, String extension) {
    index = ''
    if (extension in ['vcf.gz', 'bcf', 'bcf.gz']) {
        if (['--write-index=tbi', '-W=tbi'].any { args.contains(it) }  && extension == 'vcf.gz') {
            index = 'tbi'
        } else if (['--write-index=tbi', '-W=tbi', '--write-index=csi', '-W=csi', '--write-index', '-W'].any { args.contains(it) }) {
            index = 'csi'
        }
    }
    return index
}
