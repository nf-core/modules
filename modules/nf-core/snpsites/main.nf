process SNPSITES {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0' :
        'quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0' }"

    input:
    path alignment

    output:
    path "*.fas"        , emit: fasta
    path "*.sites.txt"  , emit: constant_sites
    env 'CONSTANT_SITES', emit: constant_sites_string
    tuple val("${task.process}"), val('snpsites'), eval('snp-sites -V 2>&1 | sed "s/snp-sites //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    snp-sites \\
        $alignment \\
        $args \\
        > filtered_alignment.fas

    echo \$(snp-sites -C $alignment) > constant.sites.txt

    export CONSTANT_SITES=\$(cat constant.sites.txt)
    """

    stub:
    """
    touch filtered_alignment.fas
    touch constant.sites.txt
    export CONSTANT_SITES=\$(cat constant.sites.txt)
    """

}
