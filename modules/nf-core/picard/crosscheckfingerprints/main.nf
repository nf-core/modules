process PICARD_CROSSCHECKFINGERPRINTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.0.0--hdfd78af_1' :
        'biocontainers/picard:3.0.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(input1)
    path input2
    path haplotype_map

    output:
    tuple val(meta), path("*.crosscheck_metrics.txt"), emit: crosscheck_metrics
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input1_string = input1.join(" --INPUT ")
    def input2_string = input2 ? "--SECOND_INPUT " + input2.join(" --SECOND_INPUT ") : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard CrosscheckFingerprints] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        CrosscheckFingerprints \\
        $args \\
        --NUM_THREADS ${task.cpus} \\
        --INPUT $input1_string \\
        $input2_string \\
        --HAPLOTYPE_MAP ${haplotype_map} \\
        --OUTPUT ${prefix}.crosscheck_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( picard CrosscheckFingerprints --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d: )
    END_VERSIONS
    """
}
