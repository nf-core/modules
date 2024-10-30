process CELLRANGERARC_MKREF {
    tag "$reference_name"
    label 'process_medium'

    container "nf-core/cellranger-arc:2.0.2"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERARC_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path fasta
    path gtf
    path motifs
    path reference_config
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "config"           , emit: config
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fast_name = fasta.name
    def gtf_name = gtf.name
    def motifs_name = motifs.name
    def reference_config = reference_config.name
    def args = task.ext.args ?: ''

    if ( !reference_name ){
        reference_name = "cellrangerarc_reference"
    }

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
    reference_config = "${reference_config}"

    if ( reference_config == "[]" ):

        config = open("config", "w")
        config.write("{\\n")
        config.write('\\torganism: "{}"\\n'.format(fasta.split(".")[0]))
        config.write('\\tgenome: ["cellrangerarc_reference"]\\n')
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangerarc: \$(echo \$( cellranger-arc --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
