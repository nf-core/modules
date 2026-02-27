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
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | sed -n "1s/.*v\\([0-9.]*\\).*/\\\\1/p"'), emit: versions_vsearch, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = user_columns ? "--userfields ${user_columns}" : ''
    if (outoption == "alnout") {
        outfmt = "--alnout"
        out_ext = 'aln'
    } else if (outoption == "biomout") {
        outfmt = "--biomout"
        out_ext = 'biom'
    } else if (outoption == "blast6out") {
        outfmt = "--blast6out"
        out_ext = 'txt'
    } else if (outoption == "mothur_shared_out") {
        outfmt = "--mothur_shared_out"
        out_ext = 'mothur'
    } else if (outoption == "otutabout") {
        outfmt = "--otutabout"
        out_ext = 'otu'
    } else if (outoption == "samout") {
        outfmt = "--samout"
        out_ext = 'sam'
    } else if (outoption == "uc") {
        outfmt = "--uc"
        out_ext = 'uc'
    } else if (outoption == "userout") {
        outfmt = "--userout"; out_ext = 'tsv'
    } else if (outoption == "lcaout") {
        outfmt = "--lcaout"
        out_ext = 'lca'
    } else {
        outfmt = "--alnout"
        out_ext = 'aln'
        log.warn("Unknown output file format provided (${outoption}): selecting pairwise alignments (alnout)")
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
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_ext = outoption == "alnout" ? 'aln' :
                    outoption == "biomout" ? 'biom' :
                    outoption == "blast6out" ? 'txt' :
                    outoption == "mothur_shared_out" ? 'mothur' :
                    outoption == "otutabout" ? 'otu' :
                    outoption == "samout" ? 'sam' :
                    outoption == "uc" ? 'uc' :
                    outoption == "userout" ? 'tsv' :
                    outoption == "lcaout" ? 'lca' :
                    'aln'
    """
    touch ${prefix}.${out_ext}
    """
}
