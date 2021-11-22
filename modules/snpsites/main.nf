process SNPSITES {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::snp-sites=2.5.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snp-sites:2.5.1--hed695b0_0"
    } else {
        container "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
    }

    input:
    path alignment

    output:
    path "*.fas"        , emit: fasta
    path "*.sites.txt"  , emit: constant_sites
    path "versions.yml" , emit: versions
    env   CONSTANT_SITES, emit: constant_sites_string

    script:
    """
    snp-sites \\
        $alignment \\
        $options.args \\
        > filtered_alignment.fas

    echo \$(snp-sites -C $alignment) > constant.sites.txt

    CONSTANT_SITES=\$(cat constant.sites.txt)

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(snp-sites -V 2>&1 | sed 's/snp-sites //')
    END_VERSIONS
    """
}
