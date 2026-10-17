process BLAST_BLASTDBCMD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/4845ad74dd459a24f390f9929801119fcf93df10161cc65292e8f416ffc05682/data'
        : 'community.wave.seqera.io/library/blast:2.17.0--9196f264d54655d4' }"

    input:
    tuple val(meta) , val(entry), path(entry_batch)
    tuple val(meta2), path(db)

    output:
    tuple val(meta), path("*.fasta"), optional: true, emit: fasta
    tuple val(meta), path("*.txt")  , optional: true, emit: text
    tuple val("${task.process}"), val('blastdbcmd'), eval('blastdbcmd -version 2>&1 | head -n1 | sed \'s/^.*blastdbcmd: //; s/ .*\$//\''), topic: versions, emit: versions_blastdbcmd

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    assert (!entry && entry_batch) || (entry && !entry_batch) : "ERROR: You must use either entry or entry_batch, not both at the same time"
    def input = ''
    if (entry) {
        input = "-entry ${entry}"
    } else {
        input = "-entry_batch ${entry_batch}"
    }
    def extension  = args.contains("-outfmt") && !args.contains("-outfmt %f") ? "txt" : "fasta"
    """
    DB=`find -L ./ -name "*.nto" | sed 's/\\.nto\$//'`
    if test -z "\$DB"
    then
        DB=`find -L ./ -name "*.pto" | sed 's/\\.pto\$//'`
    fi

    blastdbcmd \\
        -db \$DB \\
        ${args} \\
        -out ${prefix}.${extension} \\
        ${input}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension  = args.contains("-outfmt") && !args.contains("-outfmt %f") ? "txt" : "fasta"
    """
    touch ${prefix}.${extension}

    """
}
