process DSHBIO_SPLITBED {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::dsh-bio=2.0.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dsh-bio:2.0.6--hdfd78af_0' :
        'quay.io/biocontainers/dsh-bio:2.0.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    dsh-bio \\
        split-bed \\
        $args \\
        -p $prefix \\
        -s '.bed.gz' \\
        -i $bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(dsh-bio --version 2>&1 | grep -o 'dsh-bio-tools .*' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
