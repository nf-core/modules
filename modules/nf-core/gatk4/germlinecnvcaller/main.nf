process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_medium'

    if(params.enable_conda){
        error "Conda environments cannot be used for GATK4/DetermineGermlineContigPloidy at the moment. Please use docker or singularity containers."
    }
    container "broadinstitute/gatk:4.3.0.0"

    input:
    tuple val(meta), path(tsv)
    path intervals
    path model
    path ploidy

    output:
    tuple val(meta), path("*.tar.gz"), emit: tar_gz
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals_command = intervals ? "--intervals $intervals" : ""
    def model_command = model ? "--model $model" : ""
    def input_list = tsv.collect{"--input $it"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" GermlineCNVCaller \\
        $input_list \\
        --contig-ploidy-calls ${prefix}-calls \\
        --output cnv_calls/ \\
        --output-prefix $prefix \\
        $args \\
        $intervals_command \\
        $model_command
    tar -czvf cnv_calls.tar.gz cnv_calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
