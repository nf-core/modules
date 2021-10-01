// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = 0.1 // No version information printed

process CHROMAP_INDEX {
    tag '$fasta'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::chromap=0.1 bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:2cad7c5aa775241887eff8714259714a39baf016-0"
    }

    input:
    path fasta

    output:
    path "*.index"     , emit: index
    path "versions.yml", emit: version

    script:
    def prefix   = fasta.baseName
    """
    chromap \\
        -i \\
        $options.args \\
        -t $task.cpus \\
        -r $fasta \\
        -o ${prefix}.index

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo "$VERSION")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
