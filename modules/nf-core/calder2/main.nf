process CALDER2 {
    tag '$meta.id'
    label 'process_high'

    conda "bioconda::r-calder2=0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-calder2:0.3--r41hdfd78af_0' :
        'biocontainers/r-calder2:0.3--r41hdfd78af_0' }"


    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("${meta.id}/")                    , emit: output_folder
    tuple val(meta), path("${meta.id}/intermediate_data/")  , emit: intermediate_data_folder      , optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = resolution ? "::/resolutions/$resolution" : ""
    def cpus = task.cpus ?: 1
    def VERSION = '0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # getting binsize as mandatory input for calder
    binsize="\$(cooler info --field bin-size $cool$suffix)"

    calder --input $cool$suffix \\
        --outpath ${prefix} \\
        --nproc $cpus \\
        --type cool \\
        --bin_size "\${binsize}" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calder: $VERSION
    END_VERSIONS
    """
}
