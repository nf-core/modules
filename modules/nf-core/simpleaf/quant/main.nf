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

    pl_option = permitListOption(args_list, whitelist)

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
        $pl_option \\
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
// 5. 'u' (unfiltered-pl), which takes an empty string (if `chemistry` is defined as "10xv2" or "10xv3"), or a string indicating the path to a valid white list file. 
// The difference between (4) and (5) is that (4) contains the exact permit list to filter the observed barcodes, while (5) will use the white list to generate a permit list via barcode correction.

// We have two ways to take these options. `-u` is implied by the presence of the input `whitelist` channel. The options can also be passed as arguments to ext.args. Therefore, we must check two things:
// 1. if there is at least one of the options in the args list, and
// 2. if none of the four options are in the args list, there must be a non-empty whitelist channel.

def permitListOption(args_list, whitelist) {
    def pl_options = ["-k", "--knee", "-f", "--forced-cells", "-x", "--explicit-pl", "-e", "--expect-cells", "-u", "--unfiltered-pl"]
    
    // check if the args_list contains any of the pl_options
    def found = args_list.any { it in pl_options }

    // if we have an explicit pl option, we go with it and do nothing. We expect simpleaf will handle the error if there is anything wrong
    // if not, we need to check if we have an non-empty whitelist channel
    if (found) {
        return ""
    } else {
        if (whitelist) {
            return "-u <(gzip -dcf ${whitelist})"
        } else {
            error "Neither an explicit permit list generation option nor a whitelist was provided; Cannot proceed"
        }
    }
}
