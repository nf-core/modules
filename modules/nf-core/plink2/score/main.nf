process PLINK2_SCORE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::plink2=2.00a2.3" : null)
    def container_image = "plink2:2.00a2.3--h712d239_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)
    path(scorefile)

    output:
    tuple val(meta), path("*.sscore"), emit: score
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega() // plink is greedy
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --pfile ${pgen.baseName} vzs \\
        --score ${scorefile} \\
        $args \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
