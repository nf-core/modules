process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Sentieon modules does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container 'nfcore/sentieon:202112.06'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_and_bai
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [ \${SENTIEON_LICENSE_BASE64:-"unset"} != "unset" ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        if [ "\${#SENTIEON_LICENSE_BASE64}" -lt "1500" ]; then # Sentieon License server
            export SENTIEON_LICENSE=\$(echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d)
        else  # Localhost license file
            export SENTIEON_LICENSE=\$(mktemp)
            echo -e "\$LICENSE_ENCODED" | base64 -d > \$SENTIEON_LICENSE
        fi
    fi

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    sentieon bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | sentieon util sort -r $fasta -t $task.cpus -o ${prefix}.bam --sam2bam -

    touch foo.bam
    touch foo.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
