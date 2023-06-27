
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextmap=0.1.9=h9f5acd7_1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c6242a6c1a522137de7a9e9ff90779ede11cf5c5-0':
        'biocontainers/mulled-v2-f3591ce8609c7b3b33e5715333200aa5c163aa61:c6242a6c1a522137de7a9e9ff90779ede11cf5c5-0' }"

    input:
    tuple val(meta), path(input)
    path fasta

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    if [[ $input == *.pairs.gz ]]; then
        zcat $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    else
        samtools \\
            view \\
            $reference \\
            -h \\
            $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
