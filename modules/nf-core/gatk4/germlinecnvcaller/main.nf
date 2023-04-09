process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_single'

    if(params.enable_conda){
        error "Conda environments cannot be used for GATK4/GermlineCNVCaller at the moment. Please use docker or singularity containers."
    }
    container "broadinstitute/gatk:4.4.0.0"

    input:
    tuple val(meta), path(tsv)
    path intervals
    path model
    path ploidy

    output:
    tuple val(meta), path("*-calls.tar.gz"), emit: calls, optional: true
    tuple val(meta), path("*-model.tar.gz") , emit: model, optional: true
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals $intervals" : ""
    def untar_ploidy = ploidy ? (ploidy ==~ /^.*\.tar\.gz$/ ? "tar -xzf ${ploidy}" : "") : ""
    def untar_model = model ? (model ==~ /^.*\.tar\.gz$/ ? "tar -xzf ${model}" : "") : ""
    def ploidy_command = ploidy ? (ploidy ==~ /^.*\.tar\.gz$/ ? "--contig-ploidy-calls ${ploidy.toString().replace(".tar.gz","")}" : "--contig-ploidy-calls ${ploidy}") : ""
    def model_command = model ? (model ==~ /^.*\.tar\.gz$/ ? "--model ${model.toString().replace(".tar.gz","")}/${prefix}-model" : "--model ${model}/${prefix}-model") : ""
    def input_list = tsv.collect{"--input $it"}.join(' ')
    def output_command = model ? "--output ${prefix}-calls" : "--output ${prefix}-model"
    def tar_output = model ? "tar -czf ${prefix}-calls.tar.gz ${prefix}-calls" : "tar -czf ${prefix}-model.tar.gz ${prefix}-model"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
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
