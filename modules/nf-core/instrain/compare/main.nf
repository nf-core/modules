process INSTRAIN_COMPARE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::instrain=1.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/instrain:1.6.1--pyhdfd78af_0':
        'biocontainers/instrain:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(profiles)
    tuple val(meta2), path(bams)
    path stb_file

    output:
    tuple val(meta), path("*.IS_compare")   , emit: compare
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stb_args = stb_file ? "-s ${stb_file}": ''
    """
    inStrain \\
        compare \\
        -i $profiles \\
        -o ${prefix}.IS_compare \\
        --processes $task.cpus \\
        --bams $bams \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(echo \$(inStrain compare --version 2>&1) | awk 'NF{ print \$NF }')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stb_args = stb_file ? "-s ${stb_file}": ''
    """
    mkdir -p ${prefix}.IS_compare/output
    touch ${prefix}.IS_compare/output/${prefix}.IS_compare_pooled_SNV_info.tsv
    touch ${prefix}.IS_compare/output/${prefix}.IS_compare_comparisonsTable.tsv
    touch ${prefix}.IS_compare/output/${prefix}.IS_compare_pooled_SNV_data.tsv
    touch ${prefix}.IS_compare/output/${prefix}.IS_compare_pooled_SNV_data_keys.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(echo \$(inStrain compare --version 2>&1) | awk 'NF{ print \$NF }')
    END_VERSIONS
    """
}
