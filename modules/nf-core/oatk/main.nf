process OATK {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
//    conda "YOUR-TOOL-HERE"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
//        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(reads)
    path(mito_hmm)
    path(mito_hmm_h3f)
    path(mito_hmm_h3i)
    path(mito_hmm_h3m)
    path(mito_hmm_h3p)
    path(pltd_hmm)
    path(pltd_hmm_h3f)
    path(pltd_hmm_h3i)
    path(pltd_hmm_h3m)
    path(pltd_hmm_h3p)
    path(tmp)

    output:
    tuple val(meta), path("*mito.ctg.fasta")    , emit: mito_fasta, optional: true
    tuple val(meta), path("*pltd.ctg.fasta")    , emit: pltd_fasta, optional: true
    tuple val(meta), path("*mito.ctg.bed")      , emit: mito_bed, optional: true
    tuple val(meta), path("*pltd.ctg.bed")      , emit: pltd_bed, optional: true
    tuple val(meta), path("*mito.gfa")          , emit: mito_gfa, optional: true
    tuple val(meta), path("*pltd.gfa")          , emit: pltd_gfa, optional: true
    tuple val(meta), path("*annot_mito.txt")    , emit: annot_mito_txt, optional: true
    tuple val(meta), path("*annot_pltd.txt")    , emit: annot_pltd_txt, optional: true
    tuple val(meta), path("*utg.clean.gfa")     , emit: clean_gfa, optional: true
    tuple val(meta), path("*utg.final.gfa")     , emit: final_gfa, optional: true
    tuple val(meta), path("*utg.gfa")           , emit: inital_gfa, optional: true
    tuple val(meta), path("*utg.multiplex.gfa") , emit: multiplex_gfa, optional: true
    tuple val(meta), path("*utg.unzip.gfa")     , emit: unzip_gfa, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    mito_hmm_arg = ''
    if (mito_hmm) {
        mito_hmm_arg = '-m ' + mito_hmm
    }
    pltd_hmm_arg = ''
    if (pltd_hmm) {
        pltd_hmm_arg = '-p ' + pltd_hmm
    }
    tmp_arg = ''
    if (tmp) {
        tmp_arg = '-T ' + tmp
    }
    """
    oatk \\
        $args \\
        $mito_hmm_arg \\
        $pltd_hmm_arg \\
        $tmp_arg \\
        -t $task.cpus \\
        -o ${prefix} \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(oatk --version 2>&1) )
    END_VERSIONS
    """
}
