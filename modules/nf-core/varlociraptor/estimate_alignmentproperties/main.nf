process VARLOCIRAPTOR_ESTIMATE_ALIGNMENTPROPERTIES {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::varlociraptor=5.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor5.3.3--hc349b7f_0':
        'quay.io/biocontainers/varlociraptor5.3.3--hc349b7f_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    varlociraptor \\
        estimate \\
            alignment-properties \\
            $fasta \\
            $ args \\
            --bam $bam \\
            > ${prefix}_alignment-properties.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //'  ))
    END_VERSIONS
    """
}
