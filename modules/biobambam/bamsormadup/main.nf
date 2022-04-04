process BIOBAMBAM_BAMSORMADUP {
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::biobambam=2.0.183" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/biobambam:2.0.183--h9f5acd7_1' : 'quay.io/biocontainers/biobambam:2.0.183--h9f5acd7_1'}"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("*.{bam,cram}")       ,emit: bam
    tuple val(meta), path("*.bam.bai")          ,optional:true, emit: bam_index
    tuple val(meta), path("*.txt")              ,emit: metrics
    path "versions.yml"                         ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("outputformat=cram") ? "cram" : "bam"

    if (args.contains("outputformat=cram") && reference == null) error "Reference required for CRAM output."

    """
    bamsormadup \\
        $args \\
        I=$bam \\
        O=${prefix}.${suffix} \\
        M=${prefix}.txt \\
        tmpfile=$prefix \\
        threads=$task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamsormadup: \$(echo \$(bamsormadup --version 2>&1) | sed 's/^This is biobambam2 version //; s/..biobambam2 is .*\$//' )
    END_VERSIONS
    """
}
