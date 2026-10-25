process PLINK2_INDEPPAIRWISE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c7/c7dfbbefb7c4dccd8c211a1c8bf94988dcd6a58a147e5c22e2553e80af907569/data'
:         'community.wave.seqera.io/library/plink2:2.00a5.10--223badfc86e2b6e6' }"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)
    val(win)
    val(step)
    val(r2)

    output:
    tuple val(meta), path("*.prune.in")  , emit: prune_in
    tuple val(meta), path("*.prune.out") , emit: prune_out
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//'"), topic: versions, emit: versions_plink2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = plink_genotype_file.extension == 'pgen' ? '--pfile' : '--bfile'
    def input = "${plink_genotype_file.getBaseName()}"
    """
    plink2 \\
        $mode $input \\
        $args \\
        --indep-pairwise $win $step $r2 \\
        --threads $task.cpus \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prune.in
    touch ${prefix}.prune.out
    """
}
