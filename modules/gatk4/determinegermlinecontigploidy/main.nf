process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--py36hdfd78af_1' :
        'broadinstitute/gatk'}"
        //'quay.io/biocontainers/gatk4:4.2.6.1--py36hdfd78af_1'}"

    input:
    tuple val(meta), path(tsv)
    path model
    path priors

    output:
    tuple val(meta), path("ploidy.tar.gz"), emit: ploidy
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model_command = model ? "--model $model" : ""
    def priors_command = priors ? "--contig-ploidy-priors $priors" : ""
    def input_list = tsv.collect{"--input $it"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" DetermineGermlineContigPloidy \\
        $input_list \\
        --output ploidy/ \\
        --output-prefix $prefix \\
        $args \\
        $model_command \\
        $priors_command
    tar -czvf ploidy.tar.gz ploidy/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
