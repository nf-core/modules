process PYCLONEVI {
    tag "$meta.id"
    label "process_high"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyclone-vi:0.1.6--pyhdfd78af_0' :
        'biocontainers/pyclone-vi:0.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(rds_join), val(tumour_samples)

    output:
        tuple val(meta), path("*_cluster_table.csv"),   emit: ctree_input
        tuple val(meta), path("*.tsv"),                 emit: pyclone_input
        tuple val(meta), path("*_all_fits.h5"),         emit: pyclone_all_fits
        tuple val(meta), path("*_best_fit.txt"),        emit: pyclone_best_fit
        path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.py"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_cluster_table.csv
    touch ${prefix}.tsv
    touch ${prefix}_all_fits.h5
    touch ${prefix}_best_fit.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyclonevi: \$( pip show pyclone-vi | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """
}
