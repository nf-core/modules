process SATIVAEPANG_LOOPLACE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sativa-epang:0.10.0--py314hab16a5f_0' :
        'quay.io/biocontainers/sativa-epang:0.10.0--py314hab16a5f_0' }"

    input:
    tuple val(meta), path(taskdir, stageAs: "input")

    output:
    tuple val(meta), path("*.l1o_tasks"), emit: taskdir
    tuple val("${task.process}"), val('sativaepang'), eval("sativa-epang --version | cut -d' ' -f2"), topic: versions, emit: versions_sativaepang

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # loo-place writes a jplace/logs into every fold, which would mutate lootasks' own
    # output in place (nf-core/modules#12799-style -resume bug) if run directly against
    # it -- stageAs above always stages the incoming taskdir read-only as "input", so a
    # private writable copy is built here instead. fold_*/{ref.nwk,ref.fasta,query.fasta}
    # stay symlinks -- at GTDB scale these are already near-full-alignment copies per
    # fold, so duplicating them again here isn't affordable.
    mkdir "${prefix}.l1o_tasks"
    ln -s "\$(readlink -f input/manifest.json)" "${prefix}.l1o_tasks/manifest.json"
    for fold in input/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir "\$d"
        ln -s "\$(readlink -f "\$fold")"/* "\$d/"
    done

    sativa-epang \\
        -stage loo-place \\
        -taskdir "${prefix}.l1o_tasks" \\
        -T ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}.l1o_tasks"
    ln -s "\$(readlink -f input/manifest.json)" "${prefix}.l1o_tasks/manifest.json"
    for fold in input/fold_*; do
        d="${prefix}.l1o_tasks/\$(basename "\$fold")"
        mkdir "\$d"
        touch "\$d/epa_result.jplace"
    done
    """
}
