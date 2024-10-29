process CALDER2 {
    tag '$meta.id'
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-calder2:0.7--r43hdfd78af_1' :
        'biocontainers/r-calder2:0.7--r43hdfd78af_1' }"


    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("${prefix}/")                     , emit: output_folder
    tuple val(meta), path("${prefix}/intermediate_data/")   , emit: intermediate_data_folder      , optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = resolution ? "::/resolutions/$resolution" : ""
    def cpus = task.cpus ?: 1
    def VERSION = '0.7' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.7' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p ${prefix}/sub_compartments
    mkdir -p ${prefix}/sub_domains

    touch ${prefix}/sub_compartments/all_sub_compartments.bed
    touch ${prefix}/sub_compartments/all_sub_compartments.tsv
    touch ${prefix}/sub_compartments/cor_with_ref.ALL.txt
    touch ${prefix}/sub_compartments/cor_with_ref.pdf
    touch ${prefix}/sub_compartments/cor_with_ref.txt

    touch ${prefix}/sub_domains/all_nested_boundaries.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calder: $VERSION
    END_VERSIONS
    """
}
