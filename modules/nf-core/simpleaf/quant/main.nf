process SIMPLEAF_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::simpleaf=0.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.13.0--h9f5acd7_0':
        'biocontainers/simpleaf:0.13.0--h9f5acd7_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), val(chemistry), path(reads)
    path index
    val resolution
    path txp2gene
    path whitelist

    output:
    tuple val(meta), path("*_alevin_results"), emit: alevin_results
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args_list = args.tokenize()
    def prefix    = task.ext.prefix ?: "${meta.id}"

    unfiltered_command = ""
    if (whitelist) {
        unfiltered_command = "-u <(gzip -dcf ${whitelist})"
    }

    t2g_command = ""
    if (txp2gene) {
        t2g_command = "-m $txp2gene"
    }

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    simpleaf quant \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -i ${index} \\
        -c $chemistry \\
        -r $resolution \\
        -o ${prefix}_alevin_results \\
        -t $task.cpus \\
        $t2g_command \\
        $unfiltered_command \\
        $args

    [[ ! -f ${prefix}_alevin_results/af_quant/all_freq.bin ]] && cp ${prefix}_alevin_results/af_quant/permit_freq.bin ${prefix}_alevin_results/af_quant/all_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
