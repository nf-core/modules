include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KALLISTO_INDEX {
    tag '$transcript_fasta'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    tuple path (fasta)

    output:
    path tuple(meta), path(*.idx) emit: index
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)

    }
    """
    kallisto \\
        index \\
        -i ${fasta.naseName}.idx \\
        --threads $task.cpus \\
        $options.args
        $fasta
    kallisto --version | sed -e "s/kallisto //g" > ${software}.version.txt
    """
}
