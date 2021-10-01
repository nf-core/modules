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
    tuple val(meta), path(reads)
    path  hmm

    output:
    tuple val(meta), path('*.scaffolds.fa')    , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa')      , optional:true, emit: contigs
    tuple val(meta), path('*.transcripts.fa')  , optional:true, emit: transcripts
    tuple val(meta), path('*.gene_clusters.fa'), optional:true, emit: gene_clusters
    tuple val(meta), path('*.assembly.gfa')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')             , emit: log
    path  "versions.yml"                       , emit: versions

    script:
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def custom_hmms = params.spades_hmm ? "--custom-hmms $hmm" : ""
    """
    spades.py \\
        $options.args \\
        --threads $task.cpus \\
        $custom_hmms \\
        $input_reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${prefix}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${prefix}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
    fi

    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
