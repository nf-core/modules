process SOMALIER_RELATE {
    tag "somalier_relate_all"
    label 'process_low'

    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'brentp/somalier:latest' }"

    input:
    tuple val(meta), path(extract)
    path(sample_groups)
    path(ped)

    output:
    tuple val(meta), path("*somalier.res.html"),          emit: html
    tuple val(meta), path("*somalier.res.pairs.tsv"),     emit: pairsTSV
    tuple val(meta), path("*somalier.res.samples.tsv"),   emit: samplesTSV
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_list = extract.collect{"$it"}.join(' ')
    def prefix = task.ext.prefix ?: ""
    def sample_groups_command = sample_groups ? "-g $sample_groups" : ""
    def ped_command = ped ? "-p $ped" : ""

	"""
    somalier relate \
    -o ${prefix}somalier.res ${input_list} ${args} \
    ${sample_groups_command} \
    ${ped_command}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
	"""

}
