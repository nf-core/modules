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
    tuple val("${task.process}"), val('sativaepang'), eval("sed -n 's#.*share/sativa-epang-\\([0-9.]*\\)-.*#\\1#p' \$(command -v sativa-epang)"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # \$taskdir is staged from sativaepang/lootasks's own output -- a symlink under
    # local/shared-filesystem staging, but a real copy under stageInMode 'copy' (the
    # default for cloud storage without Fusion). loo-place writes a jplace (and logs)
    # inside every fold directory, which would mutate that other task's output in place
    # and make this task uncacheable across -resume (same class of bug as
    # nf-core/modules#12799). Give it a private, writable directory tree instead: fold
    # subdirectories are real (new) directories here, their read-only contents
    # (ref.nwk/ref.fasta/query.fasta) stay symlinks to avoid copying large alignments,
    # and manifest.json at the top level is a plain symlink since loo-place never writes
    # there. mv handles both staging modes: it moves a symlink as a symlink, and renames
    # a real directory in place, so nothing is deleted either way.
    mv "$taskdir" "${taskdir}.staged"
    real_taskdir=\$(readlink -f "${taskdir}.staged")
    mkdir "$taskdir"
    find "\$real_taskdir" -mindepth 1 -maxdepth 1 | while read -r entry; do
        name=\$(basename "\$entry")
        if [ -d "\$entry" ]; then
            mkdir "$taskdir/\$name"
            find "\$entry" -mindepth 1 -maxdepth 1 -exec ln -s {} "$taskdir/\$name"/ \\;
        else
            ln -s "\$entry" "$taskdir/\$name"
        fi
    done

    sativa-epang \\
        -stage loo-place \\
        -taskdir $taskdir \\
        -T ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mv "$taskdir" "${taskdir}.staged"
    real_taskdir=\$(readlink -f "${taskdir}.staged")
    mkdir "$taskdir"
    find "\$real_taskdir" -mindepth 1 -maxdepth 1 | while read -r entry; do
        name=\$(basename "\$entry")
        if [ -d "\$entry" ]; then
            mkdir "$taskdir/\$name"
            touch "$taskdir/\$name/epa_result.jplace"
        else
            ln -s "\$entry" "$taskdir/\$name"
        fi
    done
    """
}
