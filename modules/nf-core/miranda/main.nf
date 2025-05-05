process MIRANDA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miranda:3.3a--h779adbc_3':
        'biocontainers/miranda:3.3a--h779adbc_3' }"

    input:
    tuple val(meta), path(query)
    path(mirbase)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    miranda \\
        $mirbase \\
        $query \\
        $args \\
        -out ${prefix}.out

    echo "miRNA\tTarget\tScore\tEnergy_KcalMol\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tAln_len\tSubject_Identity\tQuery_Identity" > ${prefix}.txt
    grep -A 1 "Scores for this hit:" ${prefix}.out | sort | grep ">"  | cut -c 2- | tr ' ' '\t' >> ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miranda: \$(echo \$(miranda -v | sed -n 4p | sed 's/^.*miranda v//; s/microRNA.*\$//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miranda: \$(echo \$(miranda -v | sed -n 4p | sed 's/^.*miranda v//; s/microRNA.*\$//' ))
    END_VERSIONS
    """
}
