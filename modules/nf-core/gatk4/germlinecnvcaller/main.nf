process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_single'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    container "nf-core/gatk:4.4.0.0" //Biocontainers is missing a package

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "GATK4_GERMLINECNVCALLER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(tsv), path(intervals)
    path model
    path ploidy

    output:
    tuple val(meta), path("*-cnv-calls.tar.gz"), emit: calls, optional: true
    tuple val(meta), path("*-cnv-model.tar.gz"), emit: model, optional: true
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals $intervals" : ""
    def untar_ploidy = ploidy ? (ploidy.name.endsWith(".tar.gz") ? "tar -xzf ${ploidy}" : "") : ""
    def untar_model = model ? (model.name.endsWith(".tar.gz") ? "tar -xzf ${model}" : "") : ""
    def ploidy_command = ploidy ? (ploidy.name.endsWith(".tar.gz") ? "--contig-ploidy-calls ${ploidy.toString().replace(".tar.gz","")}" : "--contig-ploidy-calls ${ploidy}") : ""
    def model_command = model ? (model.name.endsWith(".tar.gz") ? "--model ${model.toString().replace(".tar.gz","")}/${prefix}-model" : "--model ${model}/${prefix}-model") : ""
    def input_list = tsv.collect{"--input $it"}.join(' ')
    def output_command = model ? "--output ${prefix}-cnv-calls" : "--output ${prefix}-cnv-model"
    def tar_output = model ? "tar -czf ${prefix}-cnv-calls.tar.gz ${prefix}-cnv-calls" : "tar -czf ${prefix}-cnv-model.tar.gz ${prefix}-cnv-model"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    ${untar_ploidy}
    ${untar_model}

    gatk --java-options "-Xmx${avail_mem}g" GermlineCNVCaller \\
        $input_list \\
        $ploidy_command \\
        $output_command \\
        --output-prefix $prefix \\
        $args \\
        $intervals_command \\
        $model_command
    ${tar_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
