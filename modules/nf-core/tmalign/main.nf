/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Module: TMALIGN
    Structural comparison of predicted candidate structures against a known
    PETase reference. Reports TM-score (normalised by both query and reference),
    RMSD, aligned length and sequence identity for every query structure.

    A TM-score above ~0.5 indicates the two structures share the same fold,
    which is the structural evidence a plastizyme candidate needs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TMALIGN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/tmalign:20240303--hd63eeec_0'
        : 'quay.io/biocontainers/tmalign:20240303--hd63eeec_0'}"

    input:
    tuple val(meta), path(pdb_structures, stageAs: 'models/*')
    path ref_pdb, stageAs: 'reference/*'

    output:
    tuple val(meta), path("*.tmalign.tsv"), emit: results
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Nextflow resolves \t and \n inside this string, so the awk program below
    // is written without any backslash escapes of its own: fields are joined
    // with OFS and each record is terminated by print's default ORS. The one
    // newline printf needs is written \\n: a real one would leave a line at
    // column 0, which stops Nextflow stripping the script's indentation, and
    // the END_VERSIONS terminator below is then never matched.
    """
    printf 'query\treference\taligned_length\ttm_score_query\ttm_score_ref\trmsd\tseq_id\\n' \
        > ${prefix}.tmalign.tsv

    for query_pdb in ${pdb_structures}; do
        query_name=\$(basename "\${query_pdb}" .pdb)

        # TM-align 20240303 labels the two entries Structure_1 (the query) and
        # Structure_2 (the reference), and prints one TM-score line for each.
        TMalign "\${query_pdb}" ${ref_pdb} ${args} \
            | awk -v OFS="\t" -v q="\${query_name}" -v r="${ref_pdb.baseName}" '
                /^Aligned length=/ { gsub(",", " "); al = \$3; rmsd = \$5; seqid = \$NF }
                /^TM-score=/       { if (tm1 == "") tm1 = \$2; else if (tm2 == "") tm2 = \$2 }
                END {
                    if (al == "") {
                        print "ERROR: could not parse TM-align output for " q | "cat 1>&2"
                        exit 1
                    }
                    print q, r, al, tm1, tm2, rmsd, seqid
                }' >> ${prefix}.tmalign.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tmalign: \$( TMalign -v 2>&1 | grep -o 'Version [0-9]*' | head -n1 | cut -d' ' -f2 )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'query\treference\taligned_length\ttm_score_query\ttm_score_ref\trmsd\tseq_id\\n' \
        > ${prefix}.tmalign.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tmalign: 20240303
    END_VERSIONS
    """
}
