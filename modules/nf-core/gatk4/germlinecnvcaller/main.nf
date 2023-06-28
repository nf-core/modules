process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "quay.io/nf-core/gatk:4.4.0.0" //Biocontainers is missing a package

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "GATK4_GERMLINECNVCALLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(tsv), path(intervals)
    tuple val(meta2), path(model)
    tuple val(meta3), path(ploidy)

    output:
    tuple val(meta), path("*-cnv-calls/*-calls"), emit: calls, optional: true
    tuple val(meta), path("*-cnv-model/*-model"), emit: model, optional: true
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals ${intervals}"         : ""
    def ploidy_command    = ploidy    ? "--contig-ploidy-calls ${ploidy}"  : ""
    def model_command     = model     ? "--model ${model}"                 : ""
    def input_list        = tsv.collect{"--input $it"}.join(' ')
    def output_command    = model     ? "--output ${prefix}-cnv-calls"     : "--output ${prefix}-cnv-model"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" GermlineCNVCaller \\
        $input_list \\
        $ploidy_command \\
        $output_command \\
        --output-prefix $prefix \\
        $args \\
        $intervals_command \\
        $model_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}-cnv-calls/${prefix}-calls
    mkdir -p ${prefix}-cnv-model/${prefix}-model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
