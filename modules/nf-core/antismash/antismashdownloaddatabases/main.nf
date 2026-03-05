process ANTISMASH_ANTISMASHDOWNLOADDATABASES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "nf-core/antismash:8.0.1--pyhdfd78af_0"

    output:
    path "antismash_db", emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download-antismash-databases \\
        --database-dir antismash_db \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "download-antismash-databases --database-dir antismash_db ${args}"

    mkdir antismash_db
    mkdir antismash_db/as-js
    mkdir antismash_db/clusterblast
    mkdir antismash_db/clustercompare
    mkdir antismash_db/comparippson
    mkdir antismash_db/knownclusterblast
    mkdir antismash_db/mite
    mkdir antismash_db/nrps_pks
    mkdir antismash_db/pfam
    mkdir antismash_db/resfam
    mkdir antismash_db/tigrfam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(echo \$(antismash --version) | sed 's/antiSMASH //;s/-.*//g')
    END_VERSIONS
    """
}
