process MITOHIFI_MITOHIFI {
    tag "$meta.id"
    label 'process_high'

    // Docker image available at the project github repository
    container 'ghcr.io/marcelauliano/mitohifi:3.2.3'

    input:
    tuple val(meta) , path(input, arity: '1..*')
    tuple val(meta2), path(ref_fa), path(ref_gb)
    val input_mode
    val mito_code

    output:
    tuple val(meta), path("final_mitogenome.fasta")         , emit: fasta
    tuple val(meta), path("contigs_stats.tsv")              , emit: stats
    tuple val(meta), path("final_mitogenome.gb")            , emit: gb                         , optional: true
    tuple val(meta), path("final_mitogenome.gff")           , emit: gff                        , optional: true
    tuple val(meta), path("all_potential_contigs.fa")       , emit: all_potential_contigs      , optional: true
    tuple val(meta), path("contigs_annotations.png")        , emit: contigs_annotations        , optional: true
    tuple val(meta), path("contigs_circularization/")       , emit: contigs_circularization    , optional: true
    tuple val(meta), path("contigs_filtering/")             , emit: contigs_filtering          , optional: true
    tuple val(meta), path("coverage_mapping/")              , emit: coverage_mapping           , optional: true
    tuple val(meta), path("coverage_plot.png")              , emit: coverage_plot              , optional: true
    tuple val(meta), path("final_mitogenome.annotation.png"), emit: final_mitogenome_annotation, optional: true
    tuple val(meta), path("final_mitogenome_choice/")       , emit: final_mitogenome_choice    , optional: true
    tuple val(meta), path("final_mitogenome.coverage.png")  , emit: final_mitogenome_coverage  , optional: true
    tuple val(meta), path("potential_contigs/")             , emit: potential_contigs          , optional: true
    tuple val(meta), path("reads_mapping_and_assembly/")    , emit: reads_mapping_and_assembly , optional: true
    tuple val(meta), path("shared_genes.tsv")               , emit: shared_genes               , optional: true
    tuple val(meta), path("*.log")                          , emit: log
    tuple val(meta), path("*")                              , emit: all_files
    // WARN: Incorrect version information is provided by tool on CLI. Please update this string when bumping container versions.
    // old version command: \$(mitohifi.py -v | sed 's/.* //')
    tuple val("${task.process}"), val('mitohifi'), eval('echo 3.2.3'), emit: versions_mitohifi, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Error: MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
    }

    // Check input compression
    def isGzipped = input.collect { f -> file(f).getExtension() == "gz" }
    if (isGzipped.any() && input_mode == "contigs") {
        error("Error: Running Mitohifi in contigs mode requires uncompressed input!")
    }
    if (isGzipped.any() && !isGzipped.every()) {
        error("Error: MitoHiFi requires all inputs to either be uncompressed or compressed!")
    }

    // Set up the input mode argument
    def modeMap = [contigs: '-c', reads: '-r']
    if (!modeMap.containsKey(input_mode)) {
        error "Error: invalid MitoHiFi input mode: ${input_mode}. Must be either 'contigs' or 'reads'!"
    }
    def input_mode_arg = modeMap[input_mode]

    // Concatenate inputs together if more than one. We have to do this via a call to
    // cat as using process substitution doesn't work as the inputs are passed to
    // other subshells
    def temp_ext = isGzipped.every() ? "fa.gz" : "fa"
    def concatenate_command = input.size() > 1 ? "cat ${input} > temp_input.${temp_ext}" : ""
    def input_string = input.size() > 1 ? "temp_input.${temp_ext}" : "${input}"
    def cleanup = input.size() > 1 ? "rm temp_input.${temp_ext}" : ""

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${concatenate_command}

    mitohifi.py \\
        ${input_mode_arg} ${input_string} \\
        -f ${ref_fa} \\
        -g ${ref_gb} \\
        -o ${mito_code} \\
        -t $task.cpus \\
        ${args} \\
        | tee ${prefix}.log

    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch final_mitogenome.fasta
    touch final_mitogenome.gb
    touch final_mitogenome.gff
    touch contigs_stats.tsv
    touch all_potential_contigs.fa
    touch contigs_annotations.png
    touch coverage_plot.png
    touch final_mitogenome.annotation.png
    touch final_mitogenome.coverage.png
    touch shared_genes.tsv
    touch shared_genes.tsv
    touch ${prefix}.log

    mkdir contigs_circularization
    mkdir contigs_filtering
    mkdir coverage_mapping
    mkdir final_mitogenome_choice
    mkdir potential_contigs
    mkdir reads_mapping_and_assembly
    """
}
