process VSEARCH_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vsearch=2.21.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'quay.io/biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(fasta)
    val clusteroption
    val idcutoff
    val outoption
    val user_columns

    output:
    tuple val(meta), path('*.aln')                , optional: true, emit: aln
    tuple val(meta), path('*.biom')               , optional: true, emit: biom
    tuple val(meta), path('*.mothur.tsv')         , optional: true, emit: mothur
    tuple val(meta), path('*.otu.tsv')            , optional: true, emit: otu
    tuple val(meta), path('*.sam')                , optional: true, emit: sam
    tuple val(meta), path('*.out.tsv')            , optional: true, emit: out
    tuple val(meta), path('*.blast.tsv')          , optional: true, emit: blast
    tuple val(meta), path('*.uc.tsv')             , optional: true, emit: uc
    tuple val(meta), path('*.centroids.fasta')    , optional: true, emit: centroids
    tuple val(meta), path('*.clusters.fasta')     , optional: true, emit: clusters
    tuple val(meta), path('*.profile.txt')        , optional: true, emit: profile
    tuple val(meta), path('*.msa.fasta')          , optional: true, emit: msa
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = user_columns ? "--userfields ${user_columns}" : ''
    switch ( clusteroption ) {
        case "fast": clustering = "--cluster_fast"; break
        case "size": clustering = "--cluster_size"; break
        case "smallmem": clustering = "--cluster_smallmem"; break
        case "unoise": clustering = "--cluster_unoise"; break
        default:
            clustering = "--cluster_fast";
            log.warn("Unknown clustering option provided (${clusteroption}): selecting fast option (--cluster_fast)");
            break
    }
    switch ( outoption ) {
        case "alnout": outfmt = "--alnout"; out_ext = 'aln'; break
        case "biomout": outfmt = "--biomout"; out_ext = 'biom'; break
        case "blast6out": outfmt = "--blast6out"; out_ext = 'blast.tsv'; break
        case "centroids": outfmt = "--centroids"; out_ext = 'centroids.fasta'; break
        case "clusters": outfmt = "--clusters"; out_ext = 'clusters.fasta'; break
        case "mothur_shared_out": outfmt = "--mothur_shared_out"; out_ext = 'mothur.tsv'; break
        case "msaout": outfmt = "--msaout"; out_ext = 'msa.fasta'; break
        case "otutabout": outfmt = "--otutabout"; out_ext = 'otu.tsv'; break
        case "profile": outfmt = "--profile"; out_ext = 'profile.txt'; break
        case "samout": outfmt = "--samout"; out_ext = 'sam'; break
        case "uc": outfmt = "--uc"; out_ext = 'uc.tsv'; break
        case "userout": outfmt = "--userout"; out_ext = 'out.tsv'; break
        default:
            outfmt = "--centroids";
            out_ext = 'centroids.fasta';
            log.warn("Unknown output file format provided (${outoption}): selecting centroids option (--centroids)");
            break
    }
    """
    vsearch \\
        ${clustering} $fasta \\
        ${outfmt} ${prefix}.${out_ext} \\
        --id $idcutoff \\
        --threads $task.cpus \\
        $args \\
        ${columns}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
