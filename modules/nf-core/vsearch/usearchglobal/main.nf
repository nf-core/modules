process VSEARCH_USEARCHGLOBAL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(queryfasta)
    path db
    val idcutoff
    val outoption
    val user_columns

    output:
    tuple val(meta), path('*.aln')    , optional: true, emit: aln
    tuple val(meta), path('*.biom')   , optional: true, emit: biom
    tuple val(meta), path('*.lca')    , optional: true, emit: lca
    tuple val(meta), path('*.mothur') , optional: true, emit: mothur
    tuple val(meta), path('*.otu')    , optional: true, emit: otu
    tuple val(meta), path('*.sam')    , optional: true, emit: sam
    tuple val(meta), path('*.tsv')    , optional: true, emit: tsv
    tuple val(meta), path('*.txt')    , optional: true, emit: txt
    tuple val(meta), path('*.uc')     , optional: true, emit: uc
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = user_columns ? "--userfields ${user_columns}" : ''
    switch ( outoption ) {
        case "alnout": outfmt = "--alnout"; out_ext = 'aln'; break
        case "biomout": outfmt = "--biomout"; out_ext = 'biom'; break
        case "blast6out": outfmt = "--blast6out"; out_ext = 'txt'; break
        case "mothur_shared_out": outfmt = "--mothur_shared_out"; out_ext = 'mothur'; break
        case "otutabout": outfmt = "--otutabout"; out_ext = 'otu'; break
        case "samout": outfmt = "--samout"; out_ext = 'sam'; break
        case "uc": outfmt = "--uc"; out_ext = 'uc'; break
        case "userout": outfmt = "--userout"; out_ext = 'tsv'; break
        case "lcaout": outfmt = "--lcaout"; out_ext = 'lca'; break
        default:
            outfmt = "--alnout";
            out_ext = 'aln';
            log.warn("Unknown output file format provided (${outoption}): selecting pairwise alignments (alnout)");
            break
    }
    """
    vsearch \\
        --usearch_global $queryfasta \\
        --db $db \\
        --id $idcutoff \\
        --threads $task.cpus \\
        $args \\
        ${columns} \\
        ${outfmt} ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
