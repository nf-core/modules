process SMOOTHXG_ITERATE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::smoothxg=0.6.7' : null)

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoothxg:0.6.7--hfb1f815_0' :
        'quay.io/biocontainers/smoothxg:0.6.7--hfb1f815_0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path("*.maf") , optional: true, emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def iterations = 2
    """
    for i in \$(seq 1 ${iterations});
    do
        input_gfa=${gfa}
        if [[ \$i != 1 ]]; then
            input_gfa=smooth.\$(echo \$(( i - 1 ))).gfa
        fi
        if [[ \$i != ${iterations} ]]; then
            smoothxg \\
                --threads=$task.cpus \\
                --gfa-in=\$input_gfa \\
                --smoothed-out=smooth.\$i.gfa \\
                $args
        else
            smoothxg \\
                --threads=$task.cpus \\
                --gfa-in=\$input_gfa \\
                --smoothed-out=${prefix}.smoothxg.gfa \\
                $args
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoothxg: \$(smoothxg --version 2>&1 | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
