process CNVKIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cnvkit=0.9.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cnvkit:0.9.9--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(tumourbam), path(normalbam)
    path  fasta
    path  targetfile

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.cnn"), emit: cnn
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*.cns"), emit: cns
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cnvkit.py \\
        batch \\
        $tumourbam \\
        --normal $normalbam\\
        --fasta $fasta \\
        --targets $targetfile \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
