process BLAST_BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/5222a42b366a0468a4c795f5057c2b8cfe39489548f8bd807e8ac0f80069bad5/data':
        'community.wave.seqera.io/library/blast:2.16.0--540f4b669b0a0ddd' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    path(taxidlist)
    val(taxids)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def taxid_mode = task.ext.taxid_mode ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta

    // Define taxid_cmd based on the args parameter
    def taxid_cmd = ""
    if (taxid_mode.contains("-taxidlist")) {
        taxid_cmd = "-taxidlist ${taxidlist}"
    } else if (taxid_mode.contains("-negative_taxidlist")) {
        taxid_cmd = "-negative_taxidlist ${taxidlist}"
    } else if (taxid_mode.contains("-taxids")) {
        taxid_cmd = "-taxids ${taxids}"
    } else if (taxid_mode.contains("-negative_taxids")) {
        taxid_cmd = "-negative_taxids ${taxids}"
    }

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.nal" | sed 's/\\.nal\$//'`
    if [ -z "\$DB" ]; then
        DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    blastn \\
        -num_threads ${task.cpus} \\
        -db \$DB \\
        -query ${fasta_name} \\
        ${taxid_cmd} \\
        ${args} \\
        -out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
