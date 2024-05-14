process MALT_BUILD {

    label 'process_high'

    conda "bioconda::malt=0.61"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/malt:0.61--hdfd78af_0' :
        'biocontainers/malt:0.61--hdfd78af_0' }"

    input:
    path fastas
    path gff
    path mapping_db

    output:
    path "malt_index/"   , emit: index
    path "versions.yml"  , emit: versions
    path "malt-build.log", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def igff = gff ? "-igff ${gff}" : ""

    """
    malt-build \\
        -v \\
        --input ${fastas.join(' ')} \\
        $igff \\
        -d 'malt_index/' \\
        -t $task.cpus \\
        $args \\
        -mdb ${mapping_db}/*.db |&tee malt-build.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        malt: \$(malt-build --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
