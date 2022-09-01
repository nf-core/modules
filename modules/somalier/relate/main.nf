process SOMALIER_RELATE {
    tag "somalier_relate_all"
    label 'process_low'

    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'brentp/somalier:latest' }"

    input:
    tuple val(meta), path(extract)

    output:
    path("*somalier.res.*"),                           emit: res
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_list = extract.collect{"$it"}.join(' ')
    def prefix = task.ext.prefix ?: ""

	"""
    somalier relate \
    -o ${prefix}somalier.res ${input_list} ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        somalier: \$(echo \$(somalier 2>&1) | sed 's/^.*somalier version: //; s/Commands:.*\$//')
    END_VERSIONS
	"""

}
