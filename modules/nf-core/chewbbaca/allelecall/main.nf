process CHEWBBACA_ALLELECALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::chewbbaca=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.1--pyhdfd78af_0' :
        'biocontainers/chewbbaca:3.3.1--pyhdfh78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(scheme)

    output:
    tuple val(meta), path("*_results_statistics.tsv"),      emit: stats
    tuple val(meta), path("*_results_contigsInfo.tsv"),     emit: contigsInfo
    tuple val(meta), path("*_results_alleles.tsv"),         emit: alleles
    tuple val(meta), path("*_paralogous_counts.tsv"),       emit: paralogous_counts, optional:true
    tuple val(meta), path("*_paralogous_loci.tsv"),         emit: paralogous_loci, optional:true
    tuple val(meta), path("*_logging_info.txt"),            emit: log
    tuple val(meta), path("*_cds_coordinates.tsv"),         emit: cds_coordinates, optional:true
    tuple val(meta), path("*_invalid_cds.txt"),             emit: invalid_cds, optional:true
    tuple val(meta), path("*_loci_summary_stats.tsv"),      emit: loci_summary_stats,   optional:true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    chewie \
      AlleleCall \
      --cpu ${task.cpus} \
      --input-files ${fasta} \
      --schema-directory ${scheme} \
      --output-directory results

    mv results/${prefix}_results_statistics.tsv ${prefix}_results_statistics.tsv
    mv results/${prefix}_results_contigsInfo.tsv ${prefix}_results_contigsInfo.tsv
    mv results/${prefix}_results_alleles.tsv ${prefix}_results_alleles.tsv
    mv results/${prefix}_paralogous_counts.tsv ${prefix}_paralogous_counts.tsv
    mv results/${prefix}_paralogous_loci.tsv ${prefix}_paralogous_loci.tsv
    mv results/${prefix}_logging_info.txt ${prefix}_logging_info.txt
    mv results/${prefix}_cds_coordinates.tsv ${prefix}_cds_coordinates.tsv
    mv results/${prefix}_invalid_cds.txt ${prefix}_invalid_cds.txt
    mv results/${prefix}_loci_summary_stats.tsv ${prefix}_loci_summary_stats.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewie //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}_results_statistics.tsv
    touch ${prefix}_results_contigsInfo.tsv
    touch ${prefix}_results_alleles.tsv
    touch ${prefix}_paralogous_counts.tsv
    touch ${prefix}_paralogous_loci.tsv
    touch ${prefix}_logging_info.txt
    touch ${prefix}_cds_coordinates.tsv
    touch ${prefix}_invalid_cds.txt
    touch ${prefix}_loci_summary_stats.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(chewie --version 2>&1) | sed 's/^.*chewie //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
