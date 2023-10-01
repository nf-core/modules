process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "nf-core/gatk:4.4.0.0" //Biocontainers is missing a package

    input:
    tuple val(meta), path(calls), path(model), path(ploidy)

    output:
    tuple val(meta), path("*_genotyped_intervals.vcf.gz") , emit: intervals, optional: true
    tuple val(meta), path("*_genotyped_segments.vcf.gz")  , emit: segments, optional: true
    tuple val(meta), path("*_denoised.vcf.gz")            , emit: denoised, optional: true
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_POSTPROCESSGERMLINECNVCALLS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def calls_command  = calls   ? calls.collect{"--calls-shard-path $it"}.join(' ')  : ""
    def model_command  = model   ? model.collect{"--model-shard-path $it"}.join(' ')  : ""
    def ploidy_command = ploidy  ? "--contig-ploidy-calls ${ploidy}"                  : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}g -XX:-UsePerfData" \\
        PostprocessGermlineCNVCalls \\
        $calls_command \\
        $model_command \\
        $ploidy_command \\
        --output-genotyped-intervals ${prefix}_genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${prefix}_genotyped_segments.vcf.gz \\
        --output-denoised-copy-ratios ${prefix}_denoised.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_POSTPROCESSGERMLINECNVCALLS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
