process PVACTOOLS_INSTALLVEPPLUGIN {
    tag 'pvactools_install_vep_plugin'
    label 'process_single'

    container "docker.io/griffithlab/pvactools:5.3.1"

    output:
    path "VEP_plugins"            , emit: vep_plugins_dir
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    git clone https://github.com/Ensembl/VEP_plugins.git
    rm -rf VEP_plugins/.git
    pvacseq install_vep_plugin VEP_plugins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvactools --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir VEP_plugins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pvactools --version)
    END_VERSIONS
    """
}
