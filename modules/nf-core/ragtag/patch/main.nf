process RAGTAG_PATCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0'
        : 'biocontainers/ragtag:2.1.0--pyhb7b1952_0'}"

    input:
    tuple val(meta), path(target), path(query)

    output:
    tuple val(meta), path("*.patched.fasta"), emit: patched_fasta
    tuple val(meta), path("*.patched.agp"),   emit: patched_agp
    tuple val(meta), path("*.comps.fasta"),   emit: patch_components_fasta
    tuple val(meta), path("*.ctg.agp"),       emit: target_splits_agp
    tuple val(meta), path("*.ctg.fasta"),     emit: target_splits_fasta
    tuple val(meta), path("*.rename.agp"),    emit: qry_rename_agp
    tuple val(meta), path("*.rename.fasta"),  emit: qry_rename_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ""
    """
    if [[ ${target} == *.gz ]]
    then
        zcat ${target} > target.fa
    else
        mv ${target} target.fa
    fi

    if [[ ${query} == *.gz ]]
    then
        zcat ${query} > query.fa
    else
        mv ${query} query.fa
    fi

    ragtag.py patch target.fa query.fa \\
        -o "${prefix}" \\
        -t ${task.cpus} \\
        ${args}

    mv ${prefix}/ragtag.patch.agp ${prefix}.patched.agp
    mv ${prefix}/ragtag.patch.fasta ${prefix}.patched.fasta
    mv ${prefix}/ragtag.patch.comps.fasta ${prefix}.comps.fasta
    mv ${prefix}/ragtag.patch.ctg.agp ${prefix}.ctg.agp
    mv ${prefix}/ragtag.patch.ctg.fasta ${prefix}.ctg.fasta
    mv ${prefix}/ragtag.patch.rename.agp ${prefix}.rename.agp
    mv ${prefix}/ragtag.patch.rename.fasta ${prefix}.rename.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RagTag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ""
    """
    touch ${prefix}.patched.agp
    touch ${prefix}.patched.fasta
    touch ${prefix}.comps.fasta
    touch ${prefix}.ctg.agp
    touch ${prefix}.ctg.fasta
    touch ${prefix}.rename.agp
    touch ${prefix}.rename.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RagTag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
