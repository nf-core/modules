process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fgbio=1.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/fgbio:1.4.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(taggedbam)
    val(strategy)

    output:
    tuple val(meta), path("*_umi-grouped.bam")  , emit: bam
    tuple val(meta), path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir tmp

    fgbio \\
        --tmp-dir=${PWD}/tmp \\
        GroupReadsByUmi \\
        -s $strategy \\
        ${options.args} \\
        -i $taggedbam \\
        -o ${prefix}_umi-grouped.bam \\
        -f ${prefix}_umi_histogram.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
