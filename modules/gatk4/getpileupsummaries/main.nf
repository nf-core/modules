process GATK4_GETPILEUPSUMMARIES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.3.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path variants
    path variants_tbi
    path sites

    output:
    tuple val(meta), path('*.pileups.table'), emit: table
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def sitesCommand = ''

    sitesCommand = sites ? " -L ${sites} " : " -L ${variants} "

    """
    gatk GetPileupSummaries \\
        -I $bam \\
        -V $variants \\
        $sitesCommand \\
        -O ${prefix}.pileups.table \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
