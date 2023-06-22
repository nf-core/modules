process CHROMOGRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::chromograph=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromograph:1.3.1--pyhdfd78af_1':
        'biocontainers/chromograph:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(autozyg)
    tuple val(meta), path(coverage)
    tuple val(meta), path(exome)
    tuple val(meta), path(fracsnp)
    tuple val(meta), path(ideogram)
    tuple val(meta), path(regions)
    tuple val(meta), path(sites)

    output:
    tuple val(meta), path("${prefix}"), emit: plots
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args   ?: ''
    prefix             = task.ext.prefix ?: "${meta.id}"
    def autozyg_param  = autozyg   ? "--autozyg ${autozyg}"   : ''
    def coverage_param = coverage  ? "--coverage ${coverage}" : ''
    def exome_param    = exome     ? "--exom ${exome}"        : ''
    def fracsnp_param  = fracsnp   ? "--fracsnp ${fracsnp}"   : ''
    def ideogram_param = ideogram  ? "--ideogram ${ideogram}" : ''
    def regions_param  = regions   ? "--regions ${regions}"   : ''
    def sites_param    = sites     ? "--sites ${sites}"       : ''
    """
    chromograph \\
        $args \\
        $autozyg_param \\
        $coverage_param \\
        $exome_param \\
        $fracsnp_param \\
        $ideogram_param \\
        $regions_param \\
        $sites_param \\
        --outd ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromograph: \$(echo \$(chromograph --version 2>&1) | sed 's/chromograph //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromograph: \$(echo \$(chromograph --version 2>&1) | sed 's/chromograph //' )
    END_VERSIONS
    """
}
