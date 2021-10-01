// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HIFIASM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::hifiasm=0.15.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hifiasm:0.15.4--h2e03b76_0"
    } else {
        container "quay.io/biocontainers/hifiasm:0.15.4--h2e03b76_0"
    }

    input:
    tuple val(meta), path(reads)
    path  paternal_kmer_dump
    path  maternal_kmer_dump
    val   use_parental_kmers

    output:
    tuple val(meta), path("*.r_utg.gfa")       , emit: raw_unitigs
    tuple val(meta), path("*.ec.bin")          , emit: corrected_reads
    tuple val(meta), path("*.ovlp.source.bin") , emit: source_overlaps
    tuple val(meta), path("*.ovlp.reverse.bin"), emit: reverse_overlaps
    tuple val(meta), path("*.p_utg.gfa")       , emit: processed_unitigs, optional: true
    tuple val(meta), path("*.asm.p_ctg.gfa")   , emit: primary_contigs  , optional: true
    tuple val(meta), path("*.asm.a_ctg.gfa")   , emit: alternate_contigs, optional: true
    tuple val(meta), path("*.hap1.p_ctg.gfa")  , emit: paternal_contigs , optional: true
    tuple val(meta), path("*.hap2.p_ctg.gfa")  , emit: maternal_contigs , optional: true
    path  "versions.yml"                       , emit: version

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (use_parental_kmers) {
        """
        hifiasm \\
            $options.args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            -1 $paternal_kmer_dump \\
            -2 $maternal_kmer_dump \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    } else { // Phasing with Hi-C data is not supported yet
        """
        hifiasm \\
            $options.args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            $reads

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    }
}
