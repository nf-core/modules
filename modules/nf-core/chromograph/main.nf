process CHROMOGRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::chromograph=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chromograph:1.3.1--pyhdfd78af_1':
        'biocontainers/chromograph:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(autozyg)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(exome)
    tuple val(meta4), path(fracsnp)
    tuple val(meta5), path(ideogram)
    tuple val(meta6), path(regions)
    tuple val(meta7), path(sites)

    output:
    tuple val(meta), path("${prefix}"), emit: plots
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args   ?: ''
    def autozyg_param  = autozyg         ? "--autozyg ${autozyg}"   : ''
    def coverage_param = coverage        ? "--coverage ${coverage}" : ''
    def exome_param    = exome           ? "--exom ${exome}"        : ''
    def fracsnp_param  = fracsnp         ? "--fracsnp ${fracsnp}"   : ''
    def ideogram_param = ideogram        ? "--ideogram ${ideogram}" : ''
    def regions_param  = regions         ? "--regions ${regions}"   : ''
    def sites_param    = sites           ? "--sites ${sites}"       : ''

    if (autozyg) {
        prefix         = task.ext.prefix ?: "${meta.id}"
    } else if (coverage) {
        prefix         = task.ext.prefix ?: "${meta2.id}"
    } else if (exome) {
        prefix         = task.ext.prefix ?: "${meta3.id}"
    } else if (fracsnp) {
        prefix         = task.ext.prefix ?: "${meta4.id}"
    } else if (ideogram) {
        prefix         = task.ext.prefix ?: "${meta5.id}"
    } else if (regions) {
        prefix         = task.ext.prefix ?: "${meta6.id}"
    } else {
        prefix         = task.ext.prefix ?: "${meta7.id}"
    }
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
    def args               = task.ext.args   ?: ''

    if (autozyg) {
        prefix             = task.ext.prefix ?: "${meta.id}"
    } else if (coverage) {
        prefix             = task.ext.prefix ?: "${meta2.id}"
    } else if (exome) {
        prefix             = task.ext.prefix ?: "${meta3.id}"
    } else if (fracsnp) {
        prefix             = task.ext.prefix ?: "${meta4.id}"
    } else if (ideogram) {
        prefix             = task.ext.prefix ?: "${meta5.id}"
    } else if (regions) {
        prefix             = task.ext.prefix ?: "${meta6.id}"
    } else {
        prefix             = task.ext.prefix ?: "${meta7.id}"
    }
    """
    mkdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromograph: \$(echo \$(chromograph --version 2>&1) | sed 's/chromograph //' )
    END_VERSIONS
    """
}
