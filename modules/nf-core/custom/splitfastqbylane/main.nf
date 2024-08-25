process CUSTOM_SPLITFASTQBYLANE {
    tag "$meta.id"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.split.fastq.gz"), emit: reads
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    read1 = [reads].flatten()[0]
    read2 = [reads].flatten().size() > 1 ? reads[1] : null
    template 'split_lanes_awk.sh'

    stub:
    """
    touch out.split.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
