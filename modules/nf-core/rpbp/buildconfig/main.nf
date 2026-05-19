process RPBP_BUILDCONFIG {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.5' :
        'quay.io/biocontainers/coreutils:9.5' }"

    input:
    tuple val(meta), path(fasta), path(gtf)
    val(genome_name)
    val(extra_yaml)

    output:
    tuple val(meta), path("rpbp_config.yaml"), emit: config
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | head -n1 | sed 's/^cat (GNU coreutils) //'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def name = genome_name ?: (meta.id ?: 'reference')
    def extras = extra_yaml ?: ''
    """
    cat > rpbp_config.yaml <<EOF
genome_base_path: rpbp_index
genome_name: ${name}
gtf: ${gtf}
fasta: ${fasta}
star_index: rpbp_index/star
ribosomal_index: rpbp_index/ribosomal
ribosomal_fasta: ${fasta}
${extras}
EOF
    """

    stub:
    """
    touch rpbp_config.yaml
    """
}
