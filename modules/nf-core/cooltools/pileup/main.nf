process COOLTOOLS_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cooltools=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.4--py39hbf8eff0_0':
        'quay.io/biocontainers/cooltools:0.5.4--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(cool), val(resolution)
    path anchor

    output:
    tuple val(meta), path("*.{npz,hdf5}"),val(format) , emit:npz
    path("versions.yml")                              , emit:versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = resolution ? "::/resolutions/$resolution" : ''
    def args_list = args.tokenize()
    format = args.contains("--out-format") ? args_list[args_list.findIndexOf{it=='--out-format'}] : 'NPZ'
    """
    cooltools pileup \\
        $args \\
        -p ${task.cpus} \\
        -o ${prefix}.${format.toLowerCase()} \\
        $cool$suffix \\
        $anchor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
