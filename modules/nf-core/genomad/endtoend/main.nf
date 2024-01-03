process GENOMAD_ENDTOEND {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.7.4--pyhdfd78af_0':
        'biocontainers/genomad:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    path  genomad_db

    output:
    tuple val(meta), path("*_aggregated_classification/*_aggregated_classification.tsv")    , emit: aggregated_classification   , optional: true
    tuple val(meta), path("*_annotate/*_taxonomy.tsv")                                      , emit: taxonomy
    tuple val(meta), path("*_find_proviruses/*_provirus.tsv")                               , emit: provirus
    tuple val(meta), path("*_score_calibration/*_compositions.tsv")                         , emit: compositions                , optional: true
    tuple val(meta), path("*_score_calibration/*_calibrated_aggregated_classification.tsv") , emit: calibrated_classification   , optional: true
    tuple val(meta), path("*_summary/*_plasmid.fna.gz")                                     , emit: plasmid_fasta
    tuple val(meta), path("*_summary/*_plasmid_genes.tsv")                                  , emit: plasmid_genes
    tuple val(meta), path("*_summary/*_plasmid_proteins.faa.gz")                            , emit: plasmid_proteins
    tuple val(meta), path("*_summary/*_plasmid_summary.tsv")                                , emit: plasmid_summary
    tuple val(meta), path("*_summary/*_virus.fna.gz")                                       , emit: virus_fasta
    tuple val(meta), path("*_summary/*_virus_genes.tsv")                                    , emit: virus_genes
    tuple val(meta), path("*_summary/*_virus_proteins.faa.gz")                              , emit: virus_proteins
    tuple val(meta), path("*_summary/*_virus_summary.tsv")                                  , emit: virus_summary
    path "versions.yml"                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        end-to-end \\
        $fasta \\
        ./ \\
        $genomad_db \\
        --threads $task.cpus \\
        $args

    gzip ./**/*.fna
    gzip ./**/*.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    """
    mkdir ${filename}_aggregated_classification
    touch ${filename}_aggregated_classification/${filename}_aggregated_classification.tsv
    mkdir ${filename}_annotate
    touch ${filename}_annotate/${filename}_taxonomy.tsv
    mkdir ${filename}_find_proviruses
    touch ${filename}_find_proviruses/${filename}_provirus.tsv
    mkdir ${filename}_marker_classification
    mkdir ${filename}_nn_classification
    mkdir ${filename}_score_calibration
    touch ${filename}_score_calibration/${filename}_calibrated_aggregated_classification.tsv
    touch ${filename}_score_calibration/${filename}_compositions.tsv
    mkdir ${filename}_summary
    touch ${filename}_summary/${filename}_plasmid.fna.gz
    touch ${filename}_summary/${filename}_plasmid_genes.tsv
    touch ${filename}_summary/${filename}_plasmid_proteins.faa.gz
    touch ${filename}_summary/${filename}_plasmid_summary.tsv
    touch ${filename}_summary/${filename}_virus.fna.gz
    touch ${filename}_summary/${filename}_virus_genes.tsv
    touch ${filename}_summary/${filename}_virus_proteins.faa.gz
    touch ${filename}_summary/${filename}_virus_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
