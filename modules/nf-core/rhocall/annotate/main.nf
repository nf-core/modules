process RHOCALL_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rhocall=0.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rhocall:0.5.1--py39hbf8eff0_0':
        'biocontainers/rhocall:0.5.1--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(roh)
    path bed

    output:
    tuple val(meta), path("*_rhocall.vcf"), emit: vcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def az_bed = bed ? "-b ${bed}" : ''
    """
    rhocall \\
        annotate \\
        $args \\
        $az_bed \\
        -r $roh \\
        -o ${prefix}_rhocall.vcf \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rhocall: \$(echo \$(rhocall --version 2>&1) | sed 's/rhocall, version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rhocall.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rhocall: \$(echo \$(rhocall --version 2>&1) | sed 's/rhocall, version //' )
    END_VERSIONS
    """
}
