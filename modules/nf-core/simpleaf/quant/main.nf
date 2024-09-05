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
    tuple val(meta4), path(whitelist)
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

    pl_option = permitListOption(args_list, whitelist) {

    }

    unfiltered_command = ""
    if (whitelist) {
        unfiltered_command = "-u <(gzip -dcf ${whitelist})"
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
        -i ${index} \\
        -1 ${forward.join( "," )} \\
        -2 ${reverse.join( "," )} \\
        -c $chemistry \\
        -r $resolution \\
        -o ${prefix} \\
        -t $task.cpus \\
        -m $txp2gene \\
        $unfiltered_command \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
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
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
    END_VERSIONS
    """
}

// We have mutual exclusive options for permit list generation. 
// 1. 'k' (knee), which is a flag for the knee method and any value provided will be ignored; 
// 2. 'f' (forced-cells), which takes an integer indicating the exact number of cells to recover; 
// 3. 'e' (expect-cells), which takes an integer indicating the expected number of cells to recover; 
// 4. 'x' (explicit-pl), which takes a string indicating the path to a valid permit list;
// 5. 'u' (unfiltered-pl), which takes an empty string (if `chemistry` is defined as "10xv2" or "10xv3"), or a string indicating the path to a valid white list file. The difference between (4) and (5) is that (4) contains the exact permit list to filter the observed barcodes, while (5) will use the white list to generate a permit list via barcode correction.

// We have two ways to take these options. `-u` is implied by the presence of the input `whitelist` channel. The other four options need to be passed as arguments to ext.args. Therefore, we must check two things:
// 1. if there is only one of the four options in the args list, and
// 2. if none of the four options are in the args list, there must be a non-empty whitelist channel.

def permitListOption(args_list, whitelist, chemistry) {
    // first, define five empty strings to hold the values of the five options
    // Notice that only one of them could be non-empty
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

    // we loop over the args list, if we see the desired flags, we record them. If we see desired flags with values, we add i by 1 and take the value.
    def i = 0
    while (i < args_list.size()) {
        // we get the flag
        def arg = args_list[i]

        if (arg == "-k" || arg == "--knee") {
            // knee doesn't take a value
            knee = '-k'
        } else if (arg == "-u" || arg == "--unfiltered-pl") {
            // we have two situations here
            // 1. if -u is the last element or it is followed by another flag(start with -), it means the chemistry must be 10xv2, 10xv3 or 10xv4 to obtain a built-in whitelist
            // 2. if -u is followed by a value, we assume that is a path to a whitelist file
            if (i == args_list.size() - 1 || args_list[i + 1].startsWith("-")) {
                // if whitelist is not provided, we check the chemistry
                // else, we use the provided whitelist
                if (whitelist) {
                    unfiltered_pl = "-u <(gzip -dcf ${whitelist})"
                } else {
                    if (chemistry == "10xv3" || chemistry == "10xv2") {
                        unfiltered_pl = '-u'
                    } else {
                        error "Unfiltered permit list is required for chemistry ${chemistry}; Cannot proceed"
                    }
                }
            } else {
                i += 1
                // now, we assume the next element is the path to the whitelist file
                wl_file = file(args_list[i], checkIfExists: true)
                unfiltered_pl = '-u ' + "${wl_file}"
            }
        } else if (arg == "-f" || arg == "--forced-cells") {
            i += 1
            // forced-cells takes an integer
            // we make sure -f is not the last element and the next element is not a flag
            if (i == args_list.size() || args_list[i].startsWith("-")) {
                error "Forced cells must be an integer; Cannot proceed"
            }
            forced_cells = '-f ' + args_list[i]
        } else if (arg == "-e" || arg == "--expect-cells") {
            i += 1
            // expect-cells takes an integer
            // we make sure -e is not the last element and the next element is not a flag
            if (i == args_list.size() || args_list[i].startsWith("-")) {
                error "Expected cells must be an integer; Cannot proceed"
            }
            expect_cells = '-p ' + args_list[i + 1]
            i += 1
        } else if (arg == "-x" || arg == "--explicit-pl") {
            i += 1
            // explicit-pl takes a file
            // we make sure -x is not the last element and the next element is not a flag
            if (i == args_list.size() || args_list[i].startsWith("-")) {
                error "Explicit permit list must be a file; Cannot proceed"
            }
            v_file = file(args_list[i], checkIfExists: true)
            explicit_pl = '-x ' + "${v_file}"
        }
        i += 1
    }

    // only one option could be non-empty
    def num_options =  (unfiltered_pl != "" ? 1 : 0) + 
            (knee != "" ? 1 : 0) + 
            (forced_cells != "" ? 1 : 0) + 
            (explicit_pl != "" ? 1 : 0) + 
            (expect_cells != "" ? 1 : 0)
    if (num_options == 1) {
        return unfiltered_pl + knee + forced_cells + explicit_pl + expect_cells
    } else {
        if (num_options == 0) {
            if (whitelist) {
                return "-u <(gzip -dcf ${whitelist})"
            } else {
                error "Neither an explicit permit list nor a whitelist was provided; Cannot proceed"

            }
        }

        error "One and only one of the permit list filtering options must be used; ${num_options} options were provided; Cannot proceed"
    }
}