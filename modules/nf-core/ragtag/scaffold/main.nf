process RAGTAG_SCAFFOLD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0'
        : 'biocontainers/ragtag:2.1.0--pyhb7b1952_0'}"

    input:
    tuple val(meta), path(assembly, name: 'assembly/*')
    tuple val(meta2), path(reference, name: 'reference/*')
    tuple val(meta3), path(exclude)
    tuple val(meta4), path(skip), path(hard_skip)

    output:
    tuple val(meta), path("*.fasta"),   emit: corrected_assembly
    tuple val(meta), path("*.agp"),     emit: corrected_agp
    tuple val(meta), path("*.stats"),   emit: corrected_stats
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_exclude = exclude ? "-e ${exclude}" : ""
    def arg_skip = skip ? "-j ${skip}" : ""
    def arg_hard_skip = hard_skip ? "-J ${hard_skip}" : ""
    """
    if [[ ${assembly} == *.gz ]]
    then
        zcat ${assembly} > assembly.fa
    else
        ln -s ${assembly} assembly.fa
    fi

    if [[ ${reference} == *.gz ]]
    then
        zcat ${reference} > reference.fa
    else
        ln -s ${reference} reference.fa
    fi

    ragtag.py scaffold reference.fa assembly.fa \\
        -o "${prefix}" \\
        -t ${task.cpus} \\
        -C \\
        ${arg_exclude} \\
        ${arg_skip} \\
        ${arg_hard_skip} \\
        ${args} \\
        2> >( tee ${prefix}.stderr.log >&2 ) \\
        | tee ${prefix}.stdout.log

    mv ${prefix}/ragtag.scaffold.fasta ${prefix}.fasta
    mv ${prefix}/ragtag.scaffold.agp ${prefix}.agp
    mv ${prefix}/ragtag.scaffold.stats ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ragtag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def _args = task.ext.args ?: ''
    def _arg_exclude = exclude ? "-e ${exclude}" : ""
    def _arg_skip = skip ? "-j ${skip}" : ""
    def _arg_hard_skip = hard_skip ? "-J ${hard_skip}" : ""
    """
    touch ${prefix}.fasta
    touch ${prefix}.agp
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
        ragtag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
