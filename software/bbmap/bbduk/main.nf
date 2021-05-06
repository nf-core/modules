include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bbmap=38.90" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bbmap:38.90--he522d1c_1"
    } else {
        container "quay.io/biocontainers/bbmap:38.90--he522d1c_1"
    }

    input:
    tuple val(meta), path(reads)
    path contaminants  
    val use_contaminants 

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path '*.version.txt'                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.trim.fastq.gz" : "out1=${prefix}_1.trim.fastq.gz out2=${prefix}_2.trim.fastq.gz"
    def contaminants_fa = use_contaminants ? "ref=$contaminants" : ''

    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ //g' | sed 's/B//g' | sed 's/G/g/g') 

    bbduk.sh \\
     -Xmx\$maxmem \\
     $raw \\
     $trimmed \\
     threads=$task.cpus \\
     $options.args \\
     $contaminants_fa \\
     &> ${prefix}.bbduk.log
    echo \$(bbversion.sh) > ${software}.version.txt

    """
}
