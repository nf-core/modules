// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ULTRA_PIPELINE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ultra_bioinformatics=0.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.0.4--pyh5e36f6f_1"
    } else {
        container "quay.io/biocontainers/ultra_bioinformatics:0.0.4--pyh5e36f6f_1"
    }

    input:
    tuple val(meta), path(reads)
    path genome
    path gtf

    output:
    tuple val(meta), path("*.sam")                              , emit: sam
    tuple val(meta), path("database.db")                        , emit: database
    tuple val(meta), path("all_splice_pairs_annotations.pickle"), emit: all_splice_pairs_annotations
    tuple val(meta), path("all_splice_sites_annotations.pickle"), emit: all_splice_sites_annotations
    tuple val(meta), path("chr_to_id.pickle")                   , emit: chr_to_id
    tuple val(meta), path("exon_choordinates_to_id.pickle")     , emit: exon_choordinates_to_id
    tuple val(meta), path("flank_choordinates.pickle")          , emit: flank_choordinates
    tuple val(meta), path("gene_to_small_segments.pickle")      , emit: gene_to_small_segments
    tuple val(meta), path("id_to_chr.pickle")                   , emit: id_to_chr
    tuple val(meta), path("max_intron_chr.pickle")              , emit: max_intron_chr
    tuple val(meta), path("parts_to_segments.pickle")           , emit: parts_to_segments
    tuple val(meta), path("ref_exon_sequences.pickle")          , emit: ref_exon_sequences
    tuple val(meta), path("ref_flank_sequences.pickle")         , emit: ref_flank_sequences
    tuple val(meta), path("ref_part_sequences.pickle")          , emit: ref_part_sequences
    tuple val(meta), path("ref_segment_sequences.pickle")       , emit: ref_segment_sequences
    tuple val(meta), path("refs_id_lengths.pickle")             , emit: refs_id_lengths
    tuple val(meta), path("refs_lengths.pickle")                , emit: refs_lengths
    tuple val(meta), path("segment_id_to_choordinates.pickle")  , emit: segment_id_to_choordinates
    tuple val(meta), path("segment_to_gene.pickle")             , emit: segment_to_gene
    tuple val(meta), path("segment_to_ref.pickle")              , emit: segment_to_ref
    tuple val(meta), path("splices_to_transcripts.pickle")      , emit: splices_to_transcripts
    tuple val(meta), path("transcripts_to_splices.pickle")      , emit: transcripts_to_splices
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    uLTRA \\
        pipeline \\
        --t $task.cpus \\
        --prefix $prefix \\
        $options.args \\
        \$(pwd)/$genome \\
        \$(pwd)/$gtf \\
        \$(pwd)/$reads \\
        ./

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( uLTRA --version|sed 's/uLTRA //g' )
    END_VERSIONS
    """
}
