process COMEBIN_RUNCOMEBIN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/comebin:1.0.4--hdfd78af_0':
        'biocontainers/comebin:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam, stageAs: "bam/*")

    output:
    tuple val(meta), path("${prefix}.comebin/comebin_res/comebin_res_bins/*"), emit: bins
    tuple val(meta), path("${prefix}.comebin/comebin_res/comebin_res.tsv")   , emit: tsv
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def clean_assembly = assembly.toString() - ".gz"
    def decompress_contigs = assembly.toString() == clean_assembly ? "" : "gunzip -q -f $assembly"
    def cleanup = decompress_contigs ? "rm ${clean_assembly}" : ""
    """
    mkdir ${prefix}.comebin
    ${decompress_contigs}

    run_comebin.sh \\
        -t ${task.cpus} \\
        -a ${clean_assembly} \\
        -p bam/ \\
        -o ${prefix}.comebin \\
        $args

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """
}
