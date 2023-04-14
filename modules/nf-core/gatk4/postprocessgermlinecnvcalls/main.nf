process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "broadinstitute/gatk:4.4.0.0" //Biocontainers is missing a package

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "GATK4_POSTPROCESSGERMLINECNVCALLS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(ploidy)
    path model
    path calls

    output:
    tuple val(meta), path("*_genotyped_intervals.vcf.gz") , emit: intervals, optional: true
    tuple val(meta), path("*_genotyped_segments.vcf.gz")  , emit: segments, optional: true
    tuple val(meta), path("*_denoised.vcf.gz")            , emit: denoised, optional: true
    path  "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def untar_ploidy = ploidy ? (ploidy.name.endsWith(".tar.gz") ? "tar -xzf ${ploidy}" : "") : ""
    def untar_model = model ? (model.name.endsWith(".tar.gz") ? "tar -xzf ${model}" : "") : ""
    def untar_calls = calls ? (calls.name.endsWith(".tar.gz") ? "tar -xzf ${calls}" : "") : ""
    def ploidy_command = ploidy ? (ploidy.name.endsWith(".tar.gz") ? "--contig-ploidy-calls ${ploidy.toString().replace(".tar.gz","")}" : "--contig-ploidy-calls ${ploidy}") : ""
    def model_command = model ? (model.name.endsWith(".tar.gz") ? "--model-shard-path ${model.toString().replace(".tar.gz","")}/${prefix}-model" : "--model-shard-path ${model}/${prefix}-model") : ""
    def calls_command = calls ? (calls.name.endsWith(".tar.gz") ? "--calls-shard-path ${calls.toString().replace(".tar.gz","")}/${prefix}-calls" : "--calls-shard-path ${model}/${prefix}-calls") : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    ${untar_ploidy}
    ${untar_model}
    ${untar_calls}

    gatk --java-options "-Xmx${avail_mem}g" PostprocessGermlineCNVCalls \\
        $ploidy_command \\
        $model_command \\
        $calls_command \\
        --output-genotyped-intervals ${prefix}_genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${prefix}_genotyped_segments.vcf.gz \\
        --output-denoised-copy-ratios ${prefix}_denoised.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genotyped_intervals.vcf.gz
    touch ${prefix}_genotyped_segments.vcf.gz
    touch ${prefix}_denoised.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
