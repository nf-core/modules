process MALT_BUILD {
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/malt:0.62--hdfd78af_0'
        : 'biocontainers/malt:0.62--hdfd78af_0'}"

    input:
    path fastas, stageAs: 'fa_folder/'
    path gff
    path mapping_db
    val map_type

    output:
    path "malt_index/", emit: index
    path "versions.yml", emit: versions
    path "malt-build.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def igff = gff ? "-igff ${gff}" : ""
    """
    INDEX=`find -L . -name '*.db' -o -name '*.abin' -type f`
    echo \$INDEX
    malt-build \\
        ${args} \\
        -v \\
        --input fa_folder/ \\
        ${igff} \\
        -d 'malt_index/' \\
        -t ${task.cpus} \\
        -${map_type} \$INDEX |&tee malt-build.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch malt-build.log
    mkdir malt_index/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-run --help  2>&1 | grep -o 'version.* ' | cut -f 1 -d ',' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
