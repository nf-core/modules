process SENTIEON_HAPLOTYPER {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Sentieon modules does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    container 'nfcore/sentieon:202112.06'

    input:
    tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model) // Not sure that input_index and dragstr_model are needed for sentieon's haplotyper. (dragstr_model is used by GATK's haplotypecaller.)
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("*.unfiltered.vcf.gz")    , optional:true, emit: vcf   // added the substring ".unfiltered" in the filename of the vcf-files since without that the g.vcf.gz-files were ending up in the vcf-channel 
    tuple val(meta), path("*.unfiltered.vcf.gz.tbi"), optional:true, emit: vcf_tbi
    tuple val(meta), path("*.g.vcf.gz")             , optional:true, emit: gvcf   // these output-files has to have the extension ".vcf.gz", otherwise the subsequent GATK-MergeVCFs will fail.
    tuple val(meta), path("*.g.vcf.gz.tbi")         , optional:true, emit: gvcf_tbi
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''    // options for the driver
    def args2 = task.ext.args2 ?: ''  // options for the vcf generation
    def args3 = task.ext.args3 ?: ''  // options for the gvcf generation
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbsnp_command = dbsnp ? "-d $dbsnp" : ""
    def interval_command = intervals ? "--interval $intervals" : ""
    def sentieon_auth_mech_base64 = task.ext.sentieon_auth_mech_base64 ?: ''
    def sentieon_auth_data_base64 = task.ext.sentieon_auth_data_base64 ?: ''

    vcf_cmd = ""
    gvcf_cmd = ""
    base_cmd = '--algo Haplotyper --genotype_model multinomial ' + dbsnp_command + ' ' // TO-DO: Avoid hardcoding genotype_module; instead move it to args2 and args3?
    // Perhaps it should be possible to use dbsnp for the vcf_cmd but not for the gvcf_cmd, or the other way round?

    if (params.haplotypecaller_emit_vcf) {  // TO-DO: Find other control mechanism for nf-core version of this module. Could we use something like say task.ext.emit_vcf instead?
        vcf_cmd = base_cmd + args2 + prefix + '.unfiltered.vcf.gz'
    }

    if (params.haplotypecaller_emit_gvcf || params.joint_germline) {  // TO-DO: Find other control mechanism for nf-core version of this module. How about task.ext.emit_gvcf ?
        gvcf_cmd = base_cmd + args3 + '--emit_mode gvcf ' + prefix + '.g.vcf.gz'
    }

    """
    export SENTIEON_LICENSE=\$(echo -n "\$SENTIEON_LICENSE_BASE64" | base64 -d)

    if  [ ${sentieon_auth_mech_base64} ] && [ ${sentieon_auth_data_base64} ]; then
        # If sentieon_auth_mech_base64 and sentieon_auth_data_base64 are non-empty strings, then Sentieon is mostly likely being run with some test-license.
        export SENTIEON_AUTH_MECH=\$(echo -n "${sentieon_auth_mech_base64}" | base64 -d)
        export SENTIEON_AUTH_DATA=\$(echo -n "${sentieon_auth_data_base64}" | base64 -d)
        echo "Decoded and exported Sentieon test-license system environment variables"
    fi

    sentieon driver $args -r $fasta -t $task.cpus -i $input $interval_command $vcf_cmd $gvcf_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unfiltered.vcf.gz
    touch ${prefix}.unfiltered.vcf.gz.tbi
    touch ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
