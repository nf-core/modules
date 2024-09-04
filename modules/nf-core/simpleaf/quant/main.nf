process SIMPLEAF_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.17.2--h919a2d8_0' :
        'biocontainers/simpleaf:0.17.2--h919a2d8_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), val(chemistry), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(txp2gene)
    tuple val(meta4)
    val resolution

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args_list = args.tokenize()
    prefix    = task.ext.prefix ?: "${meta.id}"

    unfiltered_command = ""
    if (whitelist) {
        unfiltered_command = "-u <(gzip -dcf ${whitelist})"
    }

    // expected cells
    def expect_cells = meta.expected_cells ? "--expect-cells $meta.expected_cells" : ''

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    def mapper = file("$index/piscem_idx.json").exists() ? 'salmon' : 'piscem'
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    simpleaf quant \\
        -i ${index} \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -c $chemistry \\
        -r $resolution \\
        -o ${prefix} \\
        -t $task.cpus \\
        -m $txp2gene \\
        $expect_cells \\
        $unfiltered_command \\
        $use_selective_alignment \\
        $args

    [[ ! -f ${prefix}/af_quant/all_freq.bin ]] && cp ${prefix}/af_quant/permit_freq.bin ${prefix}/af_quant/all_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        ${mapper}: \$(${mapper} --version | sed -e "s/${mapper} //g")
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    def mapper = file("$index/piscem_idx.json").exists() ? 'salmon' : 'piscem'

    """
    mkdir -p ${prefix}/af_map
    mkdir -p ${prefix}/af_quant/alevin

    touch ${prefix}/af_map/map.rad
    touch ${prefix}/af_map/unmapped_bc_count.bin
    touch ${prefix}/af_quant/alevin/quants_mat_rows.txt
    touch ${prefix}/af_quant/all_freq.bin
    touch ${prefix}/af_quant/map.collated.rad
    touch ${prefix}/af_quant/permit_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        ${mapper}: \$(${mapper} --version | sed -e "s/${mapper} //g")
    END_VERSIONS
    """
}

// define function for checking permit list generation arguments:
def check_pl_filtering_args(meta4, chemistry) {
    // -u, --unfiltered-pl
    def unfiltered_pl = ""
    // -k, --knee
    def knee = ""
    // -f, --forced-cells
    def forced_cells = ""
    // -x, explicit-pl
    def explicit_pl = ""
    // -p, --expect-cells
    def expect_cells = ""

    // now check each argument and set the corresponding param if encountered
    meta4.each { k, v ->
        if (k == "k") {
            // if knee, no need to check for value
            knee = '-k'
        } else if (k == "f") {
            // if forced-cells, value must be an integer
            if (v.isInteger()) {
                forced_cells = '-f ' + v
            } else {
                error "Forced cells must be an integer"
            }
        } else if (k == "e") {
            // if expect-cells, value must be an integer
            if (v.isInteger()) {
                expect_cells = '-p ' + v
            } else {
                error "Expected cells must be an integer"
            }
        } else if (k == "x") {
            // if explicit-pl, value must be a file
            v_file = file(v, checkIfExists: true)
            explicit_pl = '-x ' + "${v_file}"
        } else if (k == "u") {
            // if unfiltered-pl, value must be a file or an empty string
            if (v != "") {
                v_file = file(v, checkIfExists: true)
                unfiltered_pl = '-u ' + "${v_file}"
            } else if (chemistry == "10xv3" || chemistry == "10xv2") {
                unfiltered_pl = '-u <(gzip -dcf ${whitelist})'
            } else {
                error "Unfiltered permit list is required for chemistry ${chemistry}"
            }
        }
    }
    // only one option could be non-empty
    if ((unfiltered_pl != "" + knee != "" + forced_cells != "" + explicit_pl != "" + expect_cells != "") == 1) {
        return unfiltered_pl + knee + forced_cells + explicit_pl + expect_cells
    } else {
        error "Only one of the permit list filtering options can be used"
    }
}