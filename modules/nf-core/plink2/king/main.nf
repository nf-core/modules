process PLINK2_KING {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0'
        : 'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0'}"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)

    output:
    tuple val(meta), path("*.kin0"), emit: kin0
    tuple val(meta), path("*.king*"), emit: king
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = plink_genotype_file.extension == 'pgen' ? '--pfile' : '--bfile'
    def input = "${plink_genotype_file.getBaseName()}"
    """
    plink2 \\
        ${mode} ${input} \\
        --make-king-table \\
        --threads ${task.cpus} \\
        ${args} \\
        --out ${prefix}

    plink2 \\
        ${mode} ${input} \\
        --make-king triangle bin \\
        --threads ${task.cpus} \\
        ${args2} \\
        --out ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.kin0
    touch ${prefix}.king.bin
    touch ${prefix}.king.id
    """
}
