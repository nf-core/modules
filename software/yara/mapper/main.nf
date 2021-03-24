// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process YARA_MAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::yara=1.0.2 bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:f794a548b8692f29264c8984ff116c2141b90d9e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-f13549097a0d1ca36f9d4f017636fb3609f6c083:f794a548b8692f29264c8984ff116c2141b90d9e-0"
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.mapped.bam"), emit: bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if(meta.single_end) {
    """
    yara_mapper $options.args -t ${task.cpus} -f bam ${index}/yara $reads | samtools view -@ ${task.cpus} -hb -F4 > ${prefix}.mapped.bam

    echo \$(yara_mapper --help  2>&1) > ${software}.version.txt
    """
    } else {
    """
    yara_mapper $options.args -t ${task.cpus} -f bam ${index}/yara ${reads[0]} ${reads[1]} > output.bam
    samtools view -@ ${task.cpus} -hF 4 -f 0x40 -b output.bam > ${prefix}_1.mapped.bam
    samtools view -@ ${task.cpus} -hF 4 -f 0x80 -b output.bam > ${prefix}_2.mapped.bam
    echo \$(yara_mapper --version  2>&1) | grep -e "yara_mapper version:" | sed 's/yara_mapper version: //g' > ${software}.version.txt
    """
    }

}
