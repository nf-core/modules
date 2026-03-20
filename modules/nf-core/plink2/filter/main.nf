process PLINK2_FILTER {
    tag "$meta.id"
    label 'process_medium'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0':
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(plink_genotype_file), path(plink_variant_file), path(plink_sample_file)

    output:
    tuple val(meta), path("*.bim"), emit: bim, optional: true
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    tuple val(meta), path("*.fam"), emit: fam, optional: true
    tuple val(meta), path("*.pgen"),emit: pgen, optional: true
    tuple val(meta), path("*.psam"),emit: psam, optional: true
    tuple val(meta), path("*.pvar"),emit: pvar, optional: true
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = plink_genotype_file.extension == 'pgen' ? '--pfile' : '--bfile'
    def outtype = plink_genotype_file.extension == "pgen" ? '--make-pgen' : '--make-bed'
    def input = "${plink_genotype_file.getBaseName()}"
    if( "${input}" == "${prefix}" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    plink2 \\
        $mode $input \\
        $args \\
        --threads $task.cpus \\
        $outtype \\
        --out $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trio = plink_genotype_file.extension == 'pgen' ? "${prefix}.pfam ${prefix}.psam ${prefix}.pvar" : "${prefix}.bed ${prefix}.bim ${prefix}.fam"
    if( "${plink_genotype_file.getBaseName()}" == "${prefix}" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${trio}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

}
