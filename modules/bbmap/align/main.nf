// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bbmap=38.92 bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-f423f44972692ad764209668b09b84ec158bffe4:9d4ba1051c1742e6dfb8fb2ba5e474da86bba324-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-f423f44972692ad764209668b09b84ec158bffe4:9d4ba1051c1742e6dfb8fb2ba5e474da86bba324-0"
    }

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    split_cpus = Math.floor(task.cpus/2)
    input = meta.single_end ? "in=${reads}" : "in=${reads[0]} in2=${reads[1]}"

    // Set the db variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containg a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( ref.isDirectory() ) {
        if ( ref ==~ /(.\/)?ref\/?/ ) {
            db = ''
        } else {
            db = "path=${ref}"
        }
    } else {
        db = "ref=${ref}"
    }

    """
    bbmap.sh \\
        $db \\
        $input \\
        out=stdout.sam \\
        $options.args \\
        threads=$split_cpus \\
        -Xmx${task.memory.toGiga()}g | samtools view -@ ${split_cpus} $options.args2 -bhS -o ${prefix}.bam

    echo \$(bbversion.sh) > ${software}.version.txt
    """
}
