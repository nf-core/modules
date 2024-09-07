process BIOBAMBAM_BAMSORMADUP {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.183--h9f5acd7_1' : 'biocontainers/biobambam:2.0.183--h9f5acd7_1'}"

    input:
    tuple val(meta) , path(bams, stageAs: "?/*")
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam")              ,optional:true, emit: bam
    tuple val(meta), path("*.bam.bai")          ,optional:true, emit: bam_index
    tuple val(meta), path("*.cram")             ,optional:true, emit: cram
    tuple val(meta), path("*.metrics.txt")      ,emit: metrics
    path "versions.yml"                         ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("outputformat=cram") ? "cram" : "bam"
    def input_string = bams instanceof List ? bams.join(" I=") : bams
    if (args.contains("outputformat=cram") && reference == null) error "Reference required for CRAM output."

    """
    bamcat \\
        I=${input_string} \\
        level=0 \\
    | bamcollate2 \\
        level=0 \\
        $args2 \\
    | bamsormadup \\
        $args \\
        M=${prefix}.metrics.txt \\
        tmpfile=$prefix \\
        threads=$task.cpus \\
        > ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamcat: \$(echo \$(bamcat --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
        bamcollate2: \$(echo \$(bamcollate2 --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
        bamsormadup: \$(echo \$(bamsormadup --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("outputformat=cram") ? "cram" : "bam"
    if (args.contains("outputformat=cram") && reference == null) error "Reference required for CRAM output."

    """
    touch ${prefix}.${suffix}
    touch ${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamcat: \$(echo \$(bamcat --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
        bamcollate2: \$(echo \$(bamcollate2 --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
        bamsormadup: \$(echo \$(bamsormadup --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
    END_VERSIONS
    """
}
