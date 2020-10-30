// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

def VERSION = '2.29'

process BEDTOOLS_FIXBEDCOORDINATES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bedtools =2.29.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0 "
    } else {
        container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    }

    input:
        tuple val(meta), path("*.sloprefseq.bed")

    output:
        tuple val(meta), path("*.sloprefseqsorted.bed"), emit: bed
        path  "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        // sorted via chromosome, then by start position
        """
         awk -F '\\t' 'length($1) <= 5 {{ print }}' ${prefix}.sloprefseq.bed |
         awk '{{ if ($2 > $3) {{ t = $2; $2 = $3; $3 = t; }} \
         else if ($2 == $3) {{ $3 += 1; }} print $0; }}' OFS='\\t' - \
         | sort-bed > ${prefix}.sloprefseqsorted.bed
         echo $VERSION > ${software}.version.txt
        """
}
