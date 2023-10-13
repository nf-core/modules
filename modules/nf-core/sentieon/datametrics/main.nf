process SENTIEON_DATAMETRICS {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container 'nf-core/sentieon:202112.06'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path('*mq_metrics.txt') , emit: mq_metrics
    tuple val(meta), path('*qd_metrics.txt') , emit: qd_metrics
    tuple val(meta), path('*gc_summary.txt') , emit: gc_summary
    tuple val(meta), path('*gc_metrics.txt') , emit: gc_metrics
    tuple val(meta), path('*aln_metrics.txt'), emit: aln_metrics
    tuple val(meta), path('*is_metrics.txt') , emit: is_metrics
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sentieon modules do not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = bam.sort().collect{"-i $it"}.join(' ')
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''
    """
    if [ "\${#SENTIEON_LICENSE_BASE64}" -lt "1500" ]; then  # If the string SENTIEON_LICENSE_BASE64 is short, then it is an encrypted url.
        export SENTIEON_LICENSE=\$(echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d)
    else  # Localhost license file
        # The license file is stored as a nextflow variable like, for instance, this:
        # nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
        export SENTIEON_LICENSE=\$(mktemp)
        echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d > \$SENTIEON_LICENSE
    fi

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    sentieon \\
        driver \\
        -t $task.cpus \\
        -r $fasta \\
        $input \\
        $args \\
        --algo GCBias --summary ${prefix}_gc_summary.txt ${prefix}_gc_metrics.txt \\
        --algo MeanQualityByCycle ${prefix}_mq_metrics.txt \\
        --algo QualDistribution ${prefix}_qd_metrics.txt \\
        --algo InsertSizeMetricAlgo ${prefix}_is_metrics.txt  \\
        --algo AlignmentStat ${prefix}_aln_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mq_metrics.txt
    touch ${prefix}_qd_metrics.txt
    touch ${prefix}_gc_summary.txt
    touch ${prefix}_gc_metrics.txt
    touch ${prefix}_aln_metrics.txt
    touch ${prefix}_is_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
