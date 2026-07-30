process TRANSRATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c3f7b24e30eab0cde846c0d31556293afc7818492d016e6770a1fef561042405/data'
:         'community.wave.seqera.io/library/transrate_blast:092d9b5f119b9733' }"

    input:
    tuple val(meta), path(assembly)
    path reference

    output:
    tuple val(meta), path("*.assemblies.csv")                                       , emit: assemblies
    tuple val(meta), path("*.contigs.csv")                                          , emit: contigs
    tuple val(meta), path("*_mqc.csv")                                              , emit: multiqc, optional: true
    tuple val("${task.process}"), val('transrate'), eval("transrate --version"), topic: versions, emit: versions_transrate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Reference-based comparative metrics (blast+, packageable). Read-based metrics
    // (--left/--right) are deliberately not supported: they hard-require the `bam-read`
    // binary, which isn't published on any conda channel (only as a prebuilt GitHub
    // release from Blahah/transrate-tools), and their aligner (snap-aligner, pinned to
    // an old dev build) segfaults building an index for at least small test genomes.
    def reference_arg = reference ? "--reference ${reference}" : ''
    """
    transrate \\
        $args \\
        --threads $task.cpus \\
        --assembly $assembly \\
        ${reference_arg} \\
        --output .

    mv assemblies.csv ${prefix}.assemblies.csv
    mv */contigs.csv ${prefix}.contigs.csv

    # transrate writes the resolved absolute path of the input assembly into the first
    # column of the report; replace it with the sample prefix so the report content is
    # reproducible across environments and task work directories.
    sed -i "2s@^[^,]*@${prefix}@" ${prefix}.assemblies.csv

    cp ${prefix}.assemblies.csv ${prefix}_mqc.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.assemblies.csv
    touch ${prefix}.contigs.csv
    """
}
