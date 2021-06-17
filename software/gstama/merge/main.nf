// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_MERGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("merged_*.bed")              , emit: bed
    tuple path("merged_*.txt"), emit: reports

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    prefix = "merged_" + prefix
    """
    for i in *.bed
    do
        echo -e "\${i}\tcapped\t1,1,1\t\${i}" >> input.tsv
    done

    tama_merge.py -f input.tsv -d merge_dup -p $prefix $options.args
    """
}
