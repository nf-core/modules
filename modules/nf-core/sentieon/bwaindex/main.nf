process SENTIEON_BWAINDEX {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    container 'nf-core/sentieon:202112.06'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(bwa), emit: index
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "bwa/${task.ext.prefix}" : "bwa/${fasta.baseName}"
    """
    mkdir bwa

    sentieon \\
        bwa index \\
        $args \\
        -p $prefix \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir bwa

    touch bwa/genome.amb
    touch bwa/genome.ann
    touch bwa/genome.bwt
    touch bwa/genome.pac
    touch bwa/genome.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
