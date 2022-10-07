process CRUMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::crumble=0.9.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/crumble:0.9.0--hb0d9459_1':
        'quay.io/biocontainers/crumble:0.9.0--hb0d9459_1' }"

    input:
    tuple val(meta), path(input)
    val(bed)

    output:
    tuple val(meta), path("*.bam"),  emit: bam,     optional: true
    tuple val(meta), path("*.cram"), emit: cram,    optional: true
    tuple val(meta), path("*.sam"),  emit: sam,     optional: true
    tuple val(meta), path("*.bed"),  emit: bed,     optional: true
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("-O sam") ? "sam" :
                    args.contains("-O bam") ? "bam" :
                    args.contains("-O cram") ? "cram" :
                    "sam"
    def bed_output = bed ? "-b ${prefix}.bed" : ""
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    def CRUMBLE_VERSION = '0.9.0' //WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    crumble \\
        $args \\
        $bed_output \\
        $input \\
        ${prefix}.${extension}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        crumble: $CRUMBLE_VERSION
    END_VERSIONS
    """
}
