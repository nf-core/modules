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
    def sentieon_encryption_key_base64 = task.ext.sentieon_encryption_key_base64 ?: ''
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_license_message_base64 = task.ext.sentieon_license_message_base64 ?: ''

    """
    export SENTIEON_LICENSE=\$(echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d)
    # If sentieon_encryption_key_base64, sentieon_auth_mech_base64 and sentieon_license_message_base64 are set, then use those to generating corresponding env variables. (Needed for using nf-cores test-license for Sentieon.)
    if [ ${sentieon_encryption_key_base64} ] && [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_license_message_base64} ]; then
        echo "sentieon_encryption_key_base64, sentieon_auth_mech_base64 and sentieon_license_message_base64 were set"
        export SENTIEON_ENCRYPTION_KEY_BASE64=\$(echo -e "${sentieon_encryption_key_base64}" | base64 -d)
        export SENTIEON_AUTH_MECH_BASE64=\$(echo -e "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_LICENSE_MESSAGE_BASE64=\$(echo -e "${sentieon_license_message_base64}" | base64 -d)
        python3 --version
        touch foo.bam
        touch foo.bam.bai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
