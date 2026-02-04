process GEMMA_LMM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gemma:0.98.5--ha36d3ea_0':
        'community.wave.seqera.io/library/gemma:0.98.5--87bf3eea4b1ea0ad' }"

    input:
    tuple val(meta), path(genotype)
    tuple val(meta2), path(phenotype)
    tuple val(meta3), path(annot)
    tuple val(meta4), path(cXX)

    output:
    tuple val(meta), path("output/${meta.id}.out.assoc.txt"), emit: matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gemma \\
        $args \\
        -g $genotype \\
        -p $phenotype  \\
        -n 4 \\
        -a $annot \\
        -k $cXX \\
        -lmm \\
        -o ${meta.id}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir output
    touch output/${meta.id}.out.assoc.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gemma: \$(gemma --version)
    END_VERSIONS
    """
}
