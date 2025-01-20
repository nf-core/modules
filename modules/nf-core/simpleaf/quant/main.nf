process SIMPLEAF_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.18.4--ha6fb395_1':
        'biocontainers/simpleaf:0.18.4--ha6fb395_1' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), val(chemistry), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(txp2gene)
    val resolution
    tuple val(meta4), path(whitelist)
    tuple val(meta5), path(map_dir)

    output:
    tuple val(meta_out), path("${prefix}"), emit: simpleaf
    tuple val(meta_out), path(map_dir), emit: map
    tuple val(meta_out), path(quant_dir), emit: quant
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args_list = args.tokenize()
    prefix    = task.ext.prefix ?: "${meta.id}"

    if ( map_dir ) {
        mapping_args = " --map-dir ${map_dir}"
        meta_out = meta5
    } else {
        def (forward, reverse) = reads.collate(2).transpose()
        mapping_args = " -i ${index} -c ${chemistry} -1 ${forward.join( "," )} -2 ${reverse.join( "," )}"
        meta_out = meta
        map_dir = "${prefix}/af_map"
    }

    // if no whitelist is provided, we hope there will be one pl option in the args list
    pl_option = permitListOption(args_list, whitelist)
    quant_dir = "${prefix}/af_quant"

    // separate forward from reverse pairs
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    simpleaf quant \\
        $mapping_args \\
        -r $resolution \\
        -o ${prefix} \\
        -t $task.cpus \\
        -m $txp2gene \\
        $pl_option \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    meta_out = []
    """
    export ALEVIN_FRY_HOME=.

    mkdir -p ${prefix}/af_map
    mkdir -p ${prefix}/af_quant/alevin

    touch ${prefix}/af_map/map.rad
    touch ${prefix}/af_map/unmapped_bc_count.bin
    touch ${prefix}/af_quant/alevin/quants_mat_rows.txt
    touch ${prefix}/af_quant/map.collated.rad
    touch ${prefix}/af_quant/permit_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
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

    // if we have a whitelist, we can use it to generate a permit list
    // otherwise, we find is an explicit permit list generation option in the args list
    //
    if (whitelist) {
        return "-u ${whitelist}" // new alevin-fry support gz whitelist file
    } else if (found) {
        //
        return ""
    } else {
        error "No permit list generation option was provided; cannot proceed."
    }
}
