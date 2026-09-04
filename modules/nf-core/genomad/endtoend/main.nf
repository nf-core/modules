process GENOMAD_ENDTOEND {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83e31b082e82e714b01040adf90f865fdd229b6d9d926525b821d207271a3922/data'
        : 'community.wave.seqera.io/library/genomad:1.12.0--27836e6e665e84b5'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(genomad_db)

    output:
    tuple val(meta), path("${prefix}/")                                                              , emit: genomad_results
    tuple val(meta), path("${prefix}/*_aggregated_classification/*_aggregated_classification.tsv")   , emit: aggregated_classification, optional: true
    tuple val(meta), path("${prefix}/*_marker_classification/*_marker_classification.tsv")           , emit: marker_classification, optional: true
    tuple val(meta), path("${prefix}/*_annotate/*_taxonomy.tsv")                                     , emit: taxonomy
    tuple val(meta), path("${prefix}/*_find_proviruses/*_provirus.tsv")                              , emit: provirus
    tuple val(meta), path("${prefix}/*_score_calibration/*_compositions.tsv")                        , emit: compositions, optional: true
    tuple val(meta), path("${prefix}/*_score_calibration/*_calibrated_aggregated_classification.tsv"), emit: calibrated_classification, optional: true
    tuple val(meta), path("${prefix}/*_summary/*_plasmid.fna.gz")                                    , emit: plasmid_fasta
    tuple val(meta), path("${prefix}/*_summary/*_plasmid_genes.tsv")                                 , emit: plasmid_genes
    tuple val(meta), path("${prefix}/*_summary/*_plasmid_proteins.faa.gz")                           , emit: plasmid_proteins
    tuple val(meta), path("${prefix}/*_summary/*_plasmid_summary.tsv")                               , emit: plasmid_summary
    tuple val(meta), path("${prefix}/*_summary/*_virus.fna.gz")                                      , emit: virus_fasta
    tuple val(meta), path("${prefix}/*_summary/*_virus_genes.tsv")                                   , emit: virus_genes
    tuple val(meta), path("${prefix}/*_summary/*_virus_proteins.faa.gz")                             , emit: virus_proteins
    tuple val(meta), path("${prefix}/*_summary/*_virus_summary.tsv")                                 , emit: virus_summary
    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/^.*geNomad, version //; s/ .*//'"), topic: versions, emit: versions_genomad
    tuple val("${task.process}"), val('genomad_db'), eval("if [ -s ${genomad_db}/version.txt ]; then cat ${genomad_db}/version.txt; else echo 'unknown'; fi"), topic: versions, emit: versions_genomad_db
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    genomad \\
        end-to-end \\
        ${fasta} \\
        ${prefix}/ \\
        ${genomad_db} \\
        --threads ${task.cpus} \\
        ${args}

    gzip ${prefix}/**/*.fna
    gzip ${prefix}/**/*.faa
    """

    stub:
    def args = task.ext.args ?: ''
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    mkdir ${prefix}
    mkdir ${prefix}/${filename}_aggregated_classification
    touch ${prefix}/${filename}_aggregated_classification/${filename}_aggregated_classification.tsv
    mkdir ${prefix}/${filename}_annotate
    touch ${prefix}/${filename}_annotate/${filename}_taxonomy.tsv
    mkdir ${prefix}/${filename}_find_proviruses
    touch ${prefix}/${filename}_find_proviruses/${filename}_provirus.tsv
    mkdir ${prefix}/${filename}_marker_classification
    touch ${prefix}/${filename}_marker_classification/${filename}_marker_classification.tsv
    mkdir ${prefix}/${filename}_nn_classification
    mkdir ${prefix}/${filename}_score_calibration
    touch ${prefix}/${filename}_score_calibration/${filename}_calibrated_aggregated_classification.tsv
    touch ${prefix}/${filename}_score_calibration/${filename}_compositions.tsv
    mkdir ${prefix}/${filename}_summary
    echo "" | gzip > ${prefix}/${filename}_summary/${filename}_plasmid.fna.gz
    touch ${prefix}/${filename}_summary/${filename}_plasmid_genes.tsv
    echo "" | gzip > ${prefix}/${filename}_summary/${filename}_plasmid_proteins.faa.gz
    touch ${prefix}/${filename}_summary/${filename}_plasmid_summary.tsv
    echo "" | gzip > ${prefix}/${filename}_summary/${filename}_virus.fna.gz
    touch ${prefix}/${filename}_summary/${filename}_virus_genes.tsv
    echo "" | gzip > ${prefix}/${filename}_summary/${filename}_virus_proteins.faa.gz
    touch ${prefix}/${filename}_summary/${filename}_virus_summary.tsv
    """
}
