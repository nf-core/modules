process METAPHLAN_MAKEDB {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.1.1--pyhdfd78af_0' :
        'biocontainers/metaphlan:4.1.1--pyhdfd78af_0' }"

    input:

    output:
    path "metaphlan_db_latest"      , emit: db
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    metaphlan \\
        --install \\
        --nproc $task.cpus \\
        --bowtie2db metaphlan_db_latest \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    mkdir metaphlan_db_latest
    touch metaphlan_db_latest/mpa_latest
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.1.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.2.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.3.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.4.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.fna.bz2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.md5
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.pkl
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.rev.1.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.rev.2.bt2
    touch metaphlan_db_latest/mpa_v31_CHOCOPhlAn_201901.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}
