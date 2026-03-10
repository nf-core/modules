process CHROMOGRAPH {
    // $meta.id can be [] because autozyg is not required, so use all ids
    tag "${[meta, meta2, meta3, meta4, meta5, meta6, meta7].collect { meta_map -> meta_map.id }.findAll().unique().join('_') ?: 'chromograph'}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/chromograph:1.3.1--pyhdfd78af_2'
        : 'biocontainers/chromograph:1.3.1--pyhdfd78af_2'}"

    input:
    tuple val(meta), path(autozyg)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(exome)
    tuple val(meta4), path(fracsnp)
    tuple val(meta5), path(ideogram)
    tuple val(meta6), path(regions)
    tuple val(meta7), path(sites)

    output:
    tuple val(meta), path("*.png"), emit: plots, optional: true
    tuple val("${task.process}"), val('chromograph'), eval("chromograph --version | sed 's/.* //'")   , topic: versions   , emit: versions_chromograph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def autozyg_param  = autozyg       ? "--autozyg ${autozyg}"   : ''
    def coverage_param = coverage      ? "--coverage ${coverage}" : ''
    def exome_param    = exome         ? "--exom ${exome}"        : ''
    def fracsnp_param  = fracsnp       ? "--fracsnp ${fracsnp}"   : ''
    def ideogram_param = ideogram      ? "--ideogram ${ideogram}" : ''
    def regions_param  = regions       ? "--regions ${regions}"   : ''
    def sites_param    = sites         ? "--sites ${sites}"       : ''

    """
    chromograph \\
        ${args} \\
        ${autozyg_param} \\
        ${coverage_param} \\
        ${exome_param} \\
        ${fracsnp_param} \\
        ${ideogram_param} \\
        ${regions_param} \\
        ${sites_param} \\
        --outd .

    """

    stub:
    def args = task.ext.args ?: ''
    euploidy = args.contains('-e') || args.contains('--euploid')

    """
    ${touchCmd(euploidy, autozyg)}
    ${touchCmd(euploidy, coverage)}
    ${touchCmd(euploidy, exome)}
    ${touchCmd(euploidy, fracsnp)}
    ${touchCmd(euploidy, ideogram)}
    ${touchCmd(euploidy, regions)}
    ${touchCmd(euploidy, sites)}

    """
}

// Helper function to generate touch commands
def touchCmd(euploidy, input_file) {
    def chrs = euploidy ? (1..22) + ['X', 'Y', 'M'] : [1]
    input_file ? chrs.collect { chr -> "touch ${input_file}_chr${chr}.png" }.join(' ') : ''
}
