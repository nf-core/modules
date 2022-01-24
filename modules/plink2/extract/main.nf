process PLINK2_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::plink2=2.00a2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a2.3--h712d239_1' :
        'quay.io/biocontainers/plink2:2.00a2.3--h712d239_1' }"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar), path(variants)

    output:
    tuple val(meta), path("*.pgen"), emit: pgen
    tuple val(meta), path("*.psam"), emit: psam
    tuple val(meta), path("*.pvar"), emit: pvar
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$pgen" == "${prefix}.pgen" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        $args \\
        --threads $task.cpus \\
        --extract $variants \\
        --make-pgen \\
        --out ${prefix}

    echo $task.memory
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
