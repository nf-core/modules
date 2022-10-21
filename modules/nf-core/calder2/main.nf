process CALDER2 {
    tag '$meta.id'
    label 'process_high'

    conda (params.enable_conda ? "bioconda::r-calder2=0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-calder2:0.3--r41hdfd78af_0' :
        'quay.io/biocontainers/r-calder2:0.3--r41hdfd78af_0' }"


    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("${meta.id}/")            , emit: folder
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def suffix = resolution ? "::$resolution" : ""
    def cpus = task.cpus ?: 1
    def VERSION = '0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # getting binsize
    binsize="\$(cooler info --field bin-size $cool$suffix)"

    echo "Detected binsize: \${binsize}"

    # calling CALDER
    calder --input $cool$suffix --outpath ${meta.id} --nproc $cpus --type cool --bin_size "\${binsize}" $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calder: $VERSION
    END_VERSIONS
    """
}
