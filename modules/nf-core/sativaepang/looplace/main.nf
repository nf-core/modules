process SATIVAEPANG_LOOPLACE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.9.3.4--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.9.3.4--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(taskdir)

    output:
    tuple val(meta), path(taskdir), emit: taskdir
    tuple val("${task.process}"), val('sativaepang'), eval("grep -m1 -oE '[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+' \$(command -v sativa-epang)"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Private writable taskdir: loo-place writes a jplace/logs into every fold, which
    # would otherwise mutate lootasks' own output in place (nf-core/modules#12799-style
    # -resume bug). fold_*/{ref.nwk,ref.fasta,query.fasta} stay symlinks -- at GTDB scale
    # these are already near-full-alignment copies per fold, so duplicating them again
    # here isn't affordable.
    mv "$taskdir" "${taskdir}.orig"
    mkdir "$taskdir"
    ln -s "\$(readlink -f "${taskdir}.orig/manifest.json")" "$taskdir/manifest.json"
    for fold in "${taskdir}.orig"/fold_*; do
        d="$taskdir/\$(basename "\$fold")"
        mkdir "\$d"
        ln -s "\$(readlink -f "\$fold")"/* "\$d/"
    done

    sativa-epang \\
        -stage loo-place \\
        -taskdir $taskdir \\
        -T ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mv "$taskdir" "${taskdir}.orig"
    mkdir "$taskdir"
    ln -s "\$(readlink -f "${taskdir}.orig/manifest.json")" "$taskdir/manifest.json"
    for fold in "${taskdir}.orig"/fold_*; do
        d="$taskdir/\$(basename "\$fold")"
        mkdir "\$d"
        touch "\$d/epa_result.jplace"
    done
    """
}
