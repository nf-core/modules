def VERSION = '0.1' // Version information not provided by tool on CLI

process CHROMAP_INDEX {
    tag '$fasta'
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::chromap=0.1 bioconda::samtools=1.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0' :
        'quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0' }"

    input:
    path fasta

    output:
    path "*.index"     , emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = fasta.baseName
    """
    chromap \\
        -i \\
        $args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromap: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
