process BCFTOOLS_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(input), path(index), path(annotations), path(annotations_index)
    path(header_lines)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.tbi")                    , emit: tbi, optional: true
    tuple val(meta), path("*.csi")                    , emit: csi, optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def header_file = header_lines ? "--header-lines ${header_lines}" : ''
    def annotations_file = annotations ? "--annotations ${annotations}" : ''
    def extension = getVcfExtension(args);
    def index_command = !index ? "bcftools index $input" : ''

    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    $index_command

    bcftools \\
        annotate \\
        $args \\
        $annotations_file \\
        $header_file \\
        --output ${prefix}.${extension} \\
        --threads $task.cpus \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = getVcfExtension(args);
    def index = getVcfIndex(args, extension);
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = index ? "touch ${prefix}.${extension}.${index}" : ""

    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
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
