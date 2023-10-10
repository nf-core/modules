process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "nf-core/gatk:4.4.0.0" //Biocontainers is missing a package

    input:
    tuple val(meta), path(tsv), path(intervals), path(ploidy), path(model)

    output:
    tuple val(meta), path("*-cnv-model/*-calls"), emit: cohortcalls, optional: true
    tuple val(meta), path("*-cnv-model/*-model"), emit: cohortmodel, optional: true
    tuple val(meta), path("*-cnv-calls/*-calls"), emit: casecalls  , optional: true
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_GERMLINECNVCALLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
    gatk --java-options "-Xmx${avail_mem}g -XX:-UsePerfData" \\
        GermlineCNVCaller \\
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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GATK4_GERMLINECNVCALLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}-cnv-calls/${prefix}-calls
    mkdir -p ${prefix}-cnv-model/${prefix}-model
    mkdir -p ${prefix}-cnv-model/${prefix}-calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
