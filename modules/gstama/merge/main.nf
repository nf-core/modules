// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GSTAMA_MERGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gs-tama:1.0.1--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/gs-tama:1.0.1--hdfd78af_0"
    }

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("merged_*.bed"), emit: bed
    path("merged_*.txt")                 , emit: reports
    path "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    prefix = "merged_" + prefix
    """
    for i in *.bed
        do
        NLINES=\$(wc -l \${i}|cut -d " " -f 1)

        if [ "\$NLINES" -gt 0 ]; then
            echo -e "\${i}\\tcapped\\t1,1,1\\t\${i}" >> input.tsv
        fi
    done

    tama_merge.py -f input.tsv -d merge_dup -p $prefix $options.args

    echo "1.0" > ${software}.version.txt
    """
}
