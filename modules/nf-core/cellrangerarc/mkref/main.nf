process CELLRANGERARC_MKREF {
    tag "$meta.id"
    label 'process_medium'

    container "nf-core/cellranger-arc:2.0.2"

    input:
    tuple val(meta), path(fasta), path(gtf), path(motifs), path(reference_config)

    output:
    tuple val(meta), path("${prefix}"), emit: reference
    tuple val(meta), path("config")   , emit: config
    tuple val("${task.process}"), val('cellrangerarc'), eval("cellranger-arc --version 2>&1 | sed 's/cellranger-arc cellranger-arc-//'"), emit: versions_cellrangerarc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def fast_name = fasta.name
    def gtf_name = gtf.name
    def motifs_name = motifs.name
    def reference_config_name = reference_config.name
    def args = task.ext.args ?: ''

    prefix = task.ext.prefix ?: "${meta.id}_reference"

    // unlike cellranger mkref and spaceranger mkref, cellranger-arc mkref is not *yet* implemented in the
    // 10x martian runtime. It is therefore not necessary to specify --localcores and --localmem
    """
    python3 <<CODE

    from os.path import exists
    import shutil

    fasta = "${fast_name}"
    gtf = "${gtf_name}"
    motifs = "${motifs_name}"
    add = "${args}"
    reference_config = "${reference_config_name}"

    if ( reference_config == "[]" ):

        config = open("config", "w")
        config.write("{\\n")
        config.write('\\torganism: "{}"\\n'.format(fasta.split(".")[0]))
        config.write('\\tgenome: ["${prefix}"]\\n')
        config.write('\\tinput_fasta: ["{}"]\\n'.format(fasta))
        config.write('\\tinput_gtf: ["{}"]\\n'.format(gtf))
        if motifs != "[]":
            config.write('\tinput_motifs: "{}"\\n'.format(motifs))
        if add != None:
            config.write(add + "\\n")
        config.write("}")
        config.close()

        print("Wrote config file")
    else:
        if ( not exists("config") ):
            shutil.move(reference_config, "config")

    CODE

    cellranger-arc \\
        mkref \\
        --config=config \\
        --nthreads=${task.cpus} \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}_reference"
    """
    mkdir -p "${prefix}/"
    touch config
    """
}
