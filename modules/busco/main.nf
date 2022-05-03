process BUSCO {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/bsuco:5.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("$meta.id/short_summary.*.txt")                                                           , emit: short_summary_txt
    tuple val(meta), path("$meta.id/short_summary.*.json")                                                          , emit: short_summary_json
    tuple val(meta), path("$meta.id/run_*/full_table.tsv")                                                          , emit: run_full_table
    tuple val(meta), path("$meta.id/run_*/short_summary.txt")                                                       , emit: run_short_summary_txt
    tuple val(meta), path("$meta.id/run_*/short_summary.json")                                                      , emit: run_short_summary_json
    tuple val(meta), path("$meta.id/run_*/missing_busco_list.tsv")                                                  , emit: run_missing_busco_list
    tuple val(meta), path("$meta.id/run_*/hmmer_output/*.out")                                       , optional:true, emit: hmmer_output
    tuple val(meta), path("$meta.id/run_*/blast_output/coordinates.tsv")                             , optional:true, emit: blast_coordinates
    tuple val(meta), path("$meta.id/run_*/blast_output/tblastn.tsv")                                 , optional:true, emit: tblastn
    tuple val(meta), path("$meta.id/run_*/blast_output/sequences/*.temp")                            , optional:true, emit: blast_sequences
    tuple val(meta), path("$meta.id/run_*/busco_sequences/single_copy_busco_sequences/*.{fna,faa}")  , optional:true, emit: single_copy_busco_sequences
    tuple val(meta), path("$meta.id/run_*/busco_sequences/multi_copy_busco_sequences/*.{fna,faa}")   , optional:true, emit: multi_copy_busco_sequences
    tuple val(meta), path("$meta.id/run_*/busco_sequences/fragmented_busco_sequences/*.{fna,faa}")   , optional:true, emit: fragmented_busco_sequences
    tuple val(meta), path("$meta.id/prodigal_output/predicted_genes/predicted.fna")                  , optional:true, emit: predicted_fna
    tuple val(meta), path("$meta.id/prodigal_output/predicted_genes/predicted.faa")                  , optional:true, emit: predicted_faa
    tuple val(meta), path("$meta.id/translated_sequences/*.faa")                                     , optional:true, emit: translated_sequences
    path "versions.yml"                                                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // handling lineage
    def lineage = meta.lineage ? "--lineage_dataset ${meta.lineage}" : ""

    """
    gzip -cdf ${fasta} > __UNCOMPRESSED_FASTA_FILE__

    export NUMEXPR_MAX_THREADS=\$((${task.cpus}*2))

    busco \\
        --in __UNCOMPRESSED_FASTA_FILE__ \\
        --mode ${meta.mode} \\
        --out $meta.id \\
        -c $task.cpus \\
        ${lineage} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
