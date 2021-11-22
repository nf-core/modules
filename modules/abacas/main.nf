process ABACAS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::abacas=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abacas:1.3.1--pl526_0"
    } else {
        container "quay.io/biocontainers/abacas:1.3.1--pl526_0"
    }

    input:
    tuple val(meta), path(scaffold)
    path  fasta

    output:
    tuple val(meta), path('*.abacas*'), emit: results
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    abacas.pl \\
        -r $fasta \\
        -q $scaffold \\
        $args \\
        -o ${prefix}.abacas

    mv nucmer.delta ${prefix}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${prefix}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${prefix}.abacas.nucmer.tiling
    mv unused_contigs.out ${prefix}.abacas.unused.contigs.out
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//')
    END_VERSIONS
    """
}
