// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPADES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::spades=3.15.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0"
    } else {
        container "quay.io/biocontainers/spades:3.15.3--h95f258a_0"
    }

    input:
    tuple val(meta), path(illumina), path(pacbio), path(nanopore)
    path  hmm

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')      , optional:true, emit: contigs
    tuple val(meta), path('*.transcripts.fa.gz')  , optional:true, emit: transcripts
    tuple val(meta), path('*.gene_clusters.fa.gz'), optional:true, emit: gene_clusters
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')                , emit: log
    path  "versions.yml"                          , emit: versions

    script:
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = params.spades_hmm ? "--custom-hmms $hmm" : ""
    """
    spades.py \\
        $options.args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $custom_hmms \\
        $illumina_reads \\
        $pacbio_reads \\
        $nanopore_reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${prefix}.scaffolds.fa
        gzip -n ${prefix}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${prefix}.contigs.fa
        gzip -n ${prefix}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
        gzip -n ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
        gzip -n ${prefix}.assembly.gfa
    fi

    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
        gzip -n ${prefix}.gene_clusters.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
