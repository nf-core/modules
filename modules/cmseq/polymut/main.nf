include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '1.0.4'

process CMSEQ_POLYMUT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::cmseq=1.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cmseq:1.0.4--pyhb7b1952_0"
    } else {
        container "quay.io/biocontainers/cmseq:1.0.4--pyhb7b1952_0"
    }

    input:
    tuple val(meta), path(bam), path(bai), path(gff), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: polymut
    path "versions.yml"                   , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fasta_refid = fasta ? "-c $fasta" : ""
    def sortindex = bai ? "" : "--sortindex"
    """
    polymut.py \\
        $options.args \\
        $sortindex \\
        $fasta_refid \\
        --gff_file $gff \\
        $bam > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo $VERSION )
    END_VERSIONS
    """
}
