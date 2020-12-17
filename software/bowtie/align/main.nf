// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0 bioconda::samtools=1.10' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:9e14e16c284d6860574cf5b624bbc44c793cb024-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:9e14e16c284d6860574cf5b624bbc44c793cb024-0'
    }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.out'), emit: log
    path  'bowtie.version.txt', emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def unaligned = params.save_unaligned ? "--un ${prefix}.unmapped" : ''
    def endedness = meta.single_end ? "$reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX=`find -L ./ -name "*.3.ebwt" | sed 's/.3.ebwt//'`
    bowtie \\
        --threads $task.cpus \\
        --sam \\
        -x \$INDEX \\
        -q \\
        $unaligned \\
        $options.args \\
        $endedness \\
        2> ${prefix}.out \\
        | samtools view $options.args2 -@ $task.cpus -bS -o ${prefix}.bam -

    bowtie --version | head -n 1 | cut -d" " -f3 > ${software}.version.txt
    """
}
