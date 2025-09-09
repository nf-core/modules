process VIREO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vireosnp:0.5.8--pyh7cba7a3_0' :
        'biocontainers/vireosnp:0.5.8--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(cell_data), val(n_donor), path(donor_file), path(vartrix_data)
    output:
    tuple val(meta), path('*_summary.tsv')        , emit: summary
    tuple val(meta), path('*_donor_ids.tsv')      , emit: donor_ids
    tuple val(meta), path('*_prob_singlet.tsv.gz'), emit: prob_singlets
    tuple val(meta), path('*_prob_doublet.tsv.gz'), emit: prob_doublets
    path 'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = cell_data       ? "-c ${cell_data}" : "--vartrixData ${vartrix_data}"

    """
    vireo \\
        $input \\
        -N ${n_donor} \\
        -d ${donor_file} \\
        -p $task.cpus \\
        -o . \\
        $args

    mv summary.tsv ${prefix}_summary.tsv
    mv donor_ids.tsv ${prefix}_donor_ids.tsv
    mv prob_singlet.tsv.gz ${prefix}_prob_singlet.tsv.gz
    mv prob_doublet.tsv.gz ${prefix}_prob_doublet.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vireo: \$(vireo | sed '1!d ; s/Welcome to vireoSNP //; s/!//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_summary.tsv
    touch ${prefix}_donor_ids.tsv
    echo "" | gzip > ${prefix}_prob_singlet.tsv.gz
    echo "" | gzip > ${prefix}_prob_doublet.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vireo: \$(vireo | sed '1!d ; s/Welcome to vireoSNP //; s/!//')
    END_VERSIONS
    """
}
