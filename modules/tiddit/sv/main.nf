process TIDDIT_SV {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::tiddit=2.12.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/tiddit:2.12.1--py38h1773678_0"
    } else {
        container "quay.io/biocontainers/tiddit:2.12.1--py38h1773678_0"
    }

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.vcf")        , emit: vcf
    tuple val(meta), path("*.ploidy.tab") , emit: ploidy
    tuple val(meta), path("*.signals.tab"), emit: signals
    path  "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def reference = fasta == "dummy_file.txt" ? "--ref $fasta" : ""
    """
    tiddit \\
        --sv \\
        $args \\
        --bam $bam \\
        $reference \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(tiddit 2>&1) | sed 's/^.*TIDDIT-//; s/ .*\$//')
    END_VERSIONS
    """
}
