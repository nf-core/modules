process MALT_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/malt:0.62--hdfd78af_0'
        : 'quay.io/biocontainers/malt:0.62--hdfd78af_0'}"

    input:
    tuple val(meta), path(fastas, stageAs: 'fa_folder/')
    path gff
    path mapping_db
    val map_type

    output:
    tuple val(meta), path("malt_index/")   , emit: index
    tuple val(meta), path("malt-build.log"), emit: log
    tuple val("${task.process}"), val("malt"), eval("malt-build --help |& sed '/version/!d;s/.*version //;s/,.*//'"), emit: versions_malt, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def igff = gff ? "-igff ${gff}" : ""
    """
    INDEX=`find -L . \\( -name '*.db' -o -name '*.abin' \\) -type f`
    echo \$INDEX
    malt-build \\
        ${args} \\
        -v \\
        --input fa_folder/ \\
        ${igff} \\
        -d 'malt_index/' \\
        -t ${task.cpus} \\
        -${map_type} \$INDEX |& tee malt-build.log
    """

    stub:
    """
    touch malt-build.log
    mkdir malt_index/
    touch malt_index/index0.idx
    touch malt_index/ref.{db,idx,inf}
    touch malt_index/table0.{db,idx}
    touch malt_index/taxonomy.{idx,map,tre}
    """
}
