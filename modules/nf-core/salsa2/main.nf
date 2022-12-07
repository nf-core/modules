process SALSA2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::salsa2=2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salsa2:2.3--py27hee3b9ab_0':
        'quay.io/biocontainers/salsa2:2.3--py27hee3b9ab_0' }"

    input:
    // Input for salsa2 : bedfile from hicpro (nf-core/hic pipeline), draft assembly and draft assembly index
    tuple val(meta), path(fasta), path(index)
    path(bed)

    output:
    tuple val(meta), path("SALSA_output/scaffolds_FINAL.fasta"), emit: fasta
    tuple val(meta), path("SALSA_output/scaffolds_FINAL.agp"), emit: agp
    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_pipeline.py \\
        $args \\
        -a $fasta \\
        -b $bed \\
        -l $index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SALSA2: $VERSION
    END_VERSIONS
    """
}
