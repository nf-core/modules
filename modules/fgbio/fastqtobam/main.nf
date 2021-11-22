process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fgbio=1.4.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/fgbio:1.4.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(reads)
    val(read_structure)

    output:
    tuple val(meta), path("*_umi_converted.bam"), emit: umibam
    path "versions.yml"                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir tmpFolder

    fgbio \\
        --tmp-dir=${PWD}/tmpFolder \\
        FastqToBam \\
        -i $reads \\
        -o "${prefix}_umi_converted.bam" \\
        --read-structures $read_structure \\
        --sample $meta.id \\
        --library $meta.id \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
