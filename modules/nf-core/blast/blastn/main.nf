process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/4845ad74dd459a24f390f9929801119fcf93df10161cc65292e8f416ffc05682/data'
        : 'community.wave.seqera.io/library/blast:2.17.0--9196f264d54655d4' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    path taxidlist
    val taxids
    val negative_tax

    output:
    tuple val(meta), path('*.txt'), emit: txt
    tuple val("${task.process}"), val("blastn"), eval("blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//'"), topic: versions, emit: versions_blastn

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def negative_tax_cmd = negative_tax ? "negative_" : ""
    def taxidlist_cmd = taxidlist ? "-${negative_tax_cmd}taxidlist ${taxidlist}" : ""
    def taxids_cmd = taxids ? "-${negative_tax_cmd}taxids ${taxids}" : ""
    if (taxidlist_cmd.any() && taxids_cmd.any()) {
        log.error("ERROR: taxidlist and taxids can not be used at the same time, choose only one argument to use for tax id filtering.")
    }

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    export BLASTDB=${db}

    DB=`find -L ./ -name "*.nal" | sed 's/\\.nal\$//'`
    if [ -z "\$DB" ]; then
        DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    blastn \\
        -num_threads ${task.cpus} \\
        -db \$DB \\
        -query ${fasta_name} \\
        ${taxidlist_cmd} \\
        ${taxids_cmd} \\
        ${args} \\
        -out ${prefix}.txt

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    """
}
