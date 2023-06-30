process STADENIOLIB_SCRAMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::staden_io_lib=1.14.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staden_io_lib:1.14.14--h0d9da7e_3' :
        'biocontainers/staden_io_lib:1.14.14--h0d9da7e_3' }"

    input:
    tuple val(meta), path(reads)
    path(fasta)
    path(fai)
    path(gzi)

    output:
    tuple val(meta), path("*.cram") ,emit: cram
    path "*.gzi"                    ,emit: gzi, optional: true
    path "versions.yml"             ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def inputformat = reads.getExtension
    def outputformat = "cram"
    if ("-O sam" in args) {
        outputformat = "sam"
    } else if ("-O bam" in args) {
        outputformat = "bam"
    }

    def reference = if fasta && fai : "--r ${fasta}" else ""
    if (outputformat == "cram" && !reference) {
        error "Cannot convert to CRAM without a reference"
    }

    def gz_index = if gzi : "--g ${gzi}" else ""
    if (outputformat == "cram" || outputformat == "sam") {
        gz_index = ""
        warning "Cannot use gzip index for CRAM or SAM output"
    }

    """
    scramble \
        $args \
        -I ${inputformat} \
        $reference \
        -t $task.cpus \
        ${reads} \
        ${prefix}.${outputformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stadeniolib: \$(echo \$(scramble -h | head -n 1 |sed 's/^.*version //'))
    END_VERSIONS
    """
}
