process SPACERANGER_MKREF {
    tag "$fasta"
    label 'process_high'

    container "nf-core/spaceranger:9c5e7dc93c32448e"

    input:
    path(fasta)
    path(gtf)
    val(reference_name)

    output:
    path("${reference_name}"), emit: reference
    tuple val("${task.process}"), val('spaceranger'), eval('spaceranger -V | sed "s/spaceranger spaceranger-//"'), emit: versions_spaceranger, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    // --localcores is passed to the martian runtime and specifies the number of allocated jobs
    // --nthreads is passed to the STAR index generation.
    // see also https://github.com/nf-core/scrnaseq/issues/329
    """
    spaceranger \\
        mkref \\
        --genome=$reference_name \\
        --fasta=$fasta \\
        --genes=$gtf \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --nthreads=${task.cpus} \\
        $args
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir -p $reference_name

    touch ${reference_name}/genome.fa
    touch ${reference_name}/genome.fa.fai
    touch ${reference_name}/reference.json
    touch ${reference_name}/Genome
    touch ${reference_name}/SA
    touch ${reference_name}/SAindex
    touch ${reference_name}/chrLength.txt
    touch ${reference_name}/chrName.txt
    touch ${reference_name}/chrNameLength.txt
    touch ${reference_name}/chrStart.txt
    touch ${reference_name}/exonGeTrInfo.tab
    touch ${reference_name}/exonInfo.tab
    touch ${reference_name}/geneInfo.tab
    touch ${reference_name}/sjdbInfo.txt
    touch ${reference_name}/sjdbList.fromGTF.out.tab
    touch ${reference_name}/sjdbList.out.tab
    touch ${reference_name}/transcriptInfo.tab

    # Create properly compressed `genes.gtf.gz` rather than just with `touch` so
    # that it works with testing without breaking; nf-test tries opening `.gz`
    # files as compressed, which doesn't work with just `touch`
    gzip -n -c /dev/null > ${reference_name}/genes.gtf.gz
    """
}
