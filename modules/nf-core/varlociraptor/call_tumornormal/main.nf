process VARLOCIRAPTOR_CALLTUMORNORMAL {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::varlociraptor=5.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor5.3.3--hc349b7f_0':
        'quay.io/biocontainers/varlociraptor5.3.3--hc349b7f_0' }"

    input:
    tuple val(meta), path(tumor)
    tuple val(meta), path(normal)


    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    varlociraptor \
        call \
            variants \
                $args \
                --output ${prefix}.bcf \
                tumor-normal \
                    $args2 \
                    --normal ${normal} \
                    --tumor ${tumor}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //'  ))
    END_VERSIONS
    """
}
