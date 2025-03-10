process STARDIST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/stardist:0.9.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.stardist.tif"), emit: mask
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args   ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"

    """
    stardist-predict2d \\
        -i $image \\
        -o . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: \$( python -m pip show --version stardist | grep "Version" | sed -e "s/Version: //g" )
        python: \$( python --version | sed -e "s/Python //g" )
        tensorflow: \$( python -m pip show --version tensorflow | grep "Version" | sed -e "s/Version: //g" )
        tifffile: \$( python -m pip show --version tifffile | grep "Version" | sed -e "s/Version: //g" )
    END_VERSIONS
    """

    stub:
    def args               = task.ext.args   ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stardist.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: \$( python -m pip show --version stardist | grep "Version" | sed -e "s/Version: //g" )
        python: \$( python --version | sed -e "s/Python //g" )
        tensorflow: \$( python -m pip show --version tensorflow | grep "Version" | sed -e "s/Version: //g" )
        tifffile: \$( python -m pip show --version tifffile | grep "Version" | sed -e "s/Version: //g" )
    END_VERSIONS
    """
}
