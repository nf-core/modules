process RAGTAG_PATCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0'
        : 'biocontainers/ragtag:2.1.0--pyhb7b1952_0'}"

    input:
    tuple val(meta), path(target, name: 'target/*') 
    tuple val(meta2), path(query, name: 'query/*')
    tuple val(meta3), path(exclude)
    tuple val(meta4), path(skip)
    
    output:
    tuple val(meta), path("*.patch.fasta"),         emit: patch_fasta
    tuple val(meta), path("*.patch.agp"),           emit: patch_agp
    tuple val(meta), path("*.comps.fasta"),         emit: patch_components_fasta
    tuple val(meta), path("*.ragtag.patch.asm.*"),  emit: assembly_alignments,      optional: true
    tuple val(meta), path("*.ctg.agp"),             emit: target_splits_agp
    tuple val(meta), path("*.ctg.fasta"),           emit: target_splits_fasta
    tuple val(meta), path("*.rename.agp"),          emit: qry_rename_agp,           optional: true
    tuple val(meta), path("*.rename.fasta"),        emit: qry_rename_fasta,         optional: true
    tuple val(meta), path("*.patch.err"),           emit: stderr
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def arg_exclude = exclude ? "-e ${exclude}" : ""
    def arg_skip = skip ? "-j ${skip}" : ""
    """
    if [[ ${target} == *.gz ]]
    then
        zcat ${target} > target.fa
    else
        ln -s ${target} target.fa
    fi

    if [[ ${query} == *.gz ]]
    then
        zcat ${query} > query.fa
    else
        ln -s ${query} query.fa
    fi

    tail -F ${prefix}/ragtag.patch.err >&2 &
    tailpid=\$!
    ragtag.py patch target.fa query.fa \\
        -o "${prefix}" \\
        -t ${task.cpus} \\
        ${arg_exclude} \\
        ${arg_skip} \\
        ${args} \\
        2> >( tee ${prefix}.stderr.log >&2 ) \\
         | tee ${prefix}.stdout.log

    kill -TERM "\$tailpid"

    mv ${prefix}/ragtag.patch.agp ${prefix}.patch.agp
    mv ${prefix}/ragtag.patch.fasta ${prefix}.patch.fasta
    mv ${prefix}/ragtag.patch.comps.fasta ${prefix}.comps.fasta
    mv ${prefix}/ragtag.patch.ctg.agp ${prefix}.ctg.agp
    mv ${prefix}/ragtag.patch.ctg.fasta ${prefix}.ctg.fasta
    if [ -f ${prefix}/ragtag.patch.rename.agp ]; then
        mv ${prefix}/ragtag.patch.rename.agp ${prefix}.rename.agp
    fi

    if [ -f ${prefix}/ragtag.patch.rename.fasta ]; then
        mv ${prefix}/ragtag.patch.rename.fasta ${prefix}.rename.fasta
    fi
    mv ${prefix}/ragtag.patch.err ${prefix}.patch.err
    # Move the assembly files from prefix folder, and add prefix
    for alignment_file in \$(ls ${prefix}/ragtag.patch.asm.*);
        do 
            mv "\$alignment_file" "\${alignment_file/${prefix}\\//${prefix}_}"
        done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ragtag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def _args = task.ext.args ?: ""
    def _arg_exclude = exclude ? "-e ${exclude}" : ""
    def _arg_skip = skip ? "-j ${skip}" : ""
    """
    touch ${prefix}.patch.agp
    touch ${prefix}.patch.fasta
    touch ${prefix}.comps.fasta
    touch ${prefix}.ctg.agp
    touch ${prefix}.ctg.fasta
    touch ${prefix}.rename.agp
    touch ${prefix}.rename.fasta
    touch ${prefix}.ragtag.patch.asm.1
    touch ${prefix}.patch.err
    
    cat <<-END_VERSIONS > versions.yml
        ragtag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
