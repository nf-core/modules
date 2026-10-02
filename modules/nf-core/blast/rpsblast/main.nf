/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Module: RPSBLAST
    Conserved-domain annotation of candidate plastizymes against NCBI CDD.

    This is the local, offline equivalent of the NCBI Batch CD-Search web
    service: CDD ships pre-formatted RPS-BLAST profile databases, so the same
    domain assignments can be computed without submitting sequences to NCBI.
    Running it locally keeps the step reproducible, cacheable by -resume and
    free of network quotas.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLAST_RPSBLAST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data'
        : 'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759'}"

    input:
    tuple val(meta), path(fasta)
    path cdd_db

    output:
    tuple val(meta), path("*.cdsearch.tsv"), emit: hits
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-evalue 0.01 -max_target_seqs 500'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def header = 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle'
    """
    # CDD_DB_DOWNLOAD records the database prefix; fall back to the full set.
    if [ -f "${cdd_db}/db_prefix.txt" ]; then
        DB_PREFIX=\$(cat ${cdd_db}/db_prefix.txt)
    else
        DB_PREFIX=Cdd
    fi

    # printf's newline is escaped twice in the module source: escaped once,
    # Nextflow would write a real newline here, leaving a line at column 0 that
    # stops it stripping the script's indentation - and the END_VERSIONS
    # terminator below would never be matched.
    printf '${header}\\n' > ${prefix}.cdsearch.tsv

    # An empty candidate FASTA is a legitimate outcome upstream, not an error:
    # emit a header-only table so the channel contract still holds.
    if [ -s ${fasta} ]; then
        rpsblast \
            -query ${fasta} \
            -db ${cdd_db}/\${DB_PREFIX} \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
            ${args} \
            >> ${prefix}.cdsearch.tsv
    else
        echo "WARN: ${fasta} is empty — no conserved-domain search performed for ${prefix}" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rpsblast: \$( rpsblast -version 2>&1 | sed -n 's/^rpsblast: //p' | sed 's/ .*//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def header = 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle'
    """
    printf '${header}\\n' > ${prefix}.cdsearch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rpsblast: 2.17.0+
    END_VERSIONS
    """
}
