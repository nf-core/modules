process VSEARCH_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0':
        'biocontainers/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.aln.gz')                , optional: true, emit: aln
    tuple val(meta), path('*.biom.gz')               , optional: true, emit: biom
    tuple val(meta), path('*.mothur.tsv.gz')         , optional: true, emit: mothur
    tuple val(meta), path('*.otu.tsv.gz')            , optional: true, emit: otu
    tuple val(meta), path('*.bam')                   , optional: true, emit: bam
    tuple val(meta), path('*.out.tsv.gz')            , optional: true, emit: out
    tuple val(meta), path('*.blast.tsv.gz')          , optional: true, emit: blast
    tuple val(meta), path('*.uc.tsv.gz')             , optional: true, emit: uc
    tuple val(meta), path('*.centroids.fasta.gz')    , optional: true, emit: centroids
    tuple val(meta), path('*.clusters.fasta*.gz')    , optional: true, emit: clusters
    tuple val(meta), path('*.profile.txt.gz')        , optional: true, emit: profile
    tuple val(meta), path('*.msa.fasta.gz')          , optional: true, emit: msa
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (!args2.contains("--cluster_fast") && !args2.contains("--cluster_size") && !args2.contains("--cluster_smallmem") && !args2.contains("--cluster_unoise") ) {
            error "Unknown clustering option provided (${args2})"
        }
    def out_ext = args3.contains("--alnout") ? "aln" :
                    args3.contains("--biomout") ? "biom" :
                    args3.contains("--blast6out") ? "blast.tsv" :
                    args3.contains("--centroids") ? "centroids.fasta" :
                    args3.contains("--clusters") ? "clusters.fasta" :
                    args3.contains("--mothur_shared_out") ? "mothur.tsv" :
                    args3.contains("--msaout") ? "msa.fasta" :
                    args3.contains("--otutabout") ? "otu.tsv" :
                    args3.contains("--profile") ? "profile.txt" :
                    args3.contains("--samout") ? "sam" :
                    args3.contains("--uc") ? "uc.tsv" :
                    args3.contains("--userout") ? "out.tsv" :
                    ""
    if (out_ext == "") { error "Unknown output file format provided (${args3})" }
    """
    vsearch \\
        $args2 $fasta \\
        $args3 ${prefix}.${out_ext} \\
        --threads $task.cpus \\
        $args

    if [[ $args3 == "--clusters" ]]
    then
        gzip -n ${prefix}.${out_ext}*
    elif [[ $args3 != "--samout" ]]
    then
        gzip -n ${prefix}.${out_ext}
    else
        samtools view -T $fasta -S -b ${prefix}.${out_ext} > ${prefix}.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
