process RHOCALL_VIZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rhocall:0.5.1--py39hbf8eff0_0':
        'biocontainers/rhocall:0.5.1--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(roh)

    output:
    tuple val(meta), path("${prefix}/${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}/${prefix}.wig"), emit: wig
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    rhocall \\
        viz \\
        $args \\
        -r $roh \\
        --out_dir ${prefix} \\
        $vcf

    mv ${prefix}/output.bed ${prefix}/${prefix}.bed
    mv ${prefix}/output.wig ${prefix}/${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rhocall: \$(echo \$(rhocall --version 2>&1) | sed 's/rhocall, version //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.bed
    touch ${prefix}/${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rhocall: \$(echo \$(rhocall --version 2>&1) | sed 's/rhocall, version //' )
    END_VERSIONS
    """
}
