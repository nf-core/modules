process VSEARCH_CLUSTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vsearch=2.21.1 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:6e91f82b987ac1c6b6a24e292a934ae41b17311d-0':
        'quay.io/biocontainers/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:6e91f82b987ac1c6b6a24e292a934ae41b17311d-0' }"

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
    tuple val(meta), path('*.clusters.fasta.gz')     , optional: true, emit: clusters
    tuple val(meta), path('*.profile.txt.gz')        , optional: true, emit: profile
    tuple val(meta), path('*.msa.fasta.gz')          , optional: true, emit: msa
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args      ?: ''
    def prefix        = task.ext.prefix    ?: "${meta.id}"

    def cluster_possible = ["--cluster_fast", "--cluster_size", "--cluster_smallmem", "--cluster_unoise"]
    def out_possible = ["--alnout", "--biomout", "--blast6out", "--centroids", "--clusters", "--mothur_shared_out", "--msaout", "--otutabout", "--profile", "--samout", "--uc", "--userout"]

    def cluster_option = args.cluster_option in cluster_possible ? "$args.cluster_option" : '--cluster_fast'
    def out_option     = args.out_option in out_possible         ? "$args.out_option"     : '--centroids'
    def id_cutoff      = args.id_cutoff                          ? "$args.id_cutoff"      : '1'
    def user_columns   = args.user_columns                       ? "$args.user_columns"   : ''
    def options        = args.options                            ? "$args.options"        : ''

    switch ( out_option ) {
        case "--alnout"            : out_ext = 'aln'; break
        case "--biomout"           : out_ext = 'biom'; break
        case "--blast6out"         : out_ext = 'blast.tsv'; break
        case "--centroids"         : out_ext = 'centroids.fasta'; break
        case "--clusters"          : out_ext = 'clusters.fasta'; break
        case "--mothur_shared_out" : out_ext = 'mothur.tsv'; break
        case "--msaout"            : out_ext = 'msa.fasta'; break
        case "--otutabout"         : out_ext = 'otu.tsv'; break
        case "--profile"           : out_ext = 'profile.txt'; break
        case "--samout"            : out_ext = 'sam'; break
        case "--uc"                : out_ext = 'uc.tsv'; break
        case "--userout"           : out_ext = 'out.tsv'; break
    }
    """
    vsearch \\
        ${cluster_option} $fasta \\
        ${out_option} ${prefix}.${out_ext} \\
        --id $id_cutoff \\
        --threads $task.cpus \\
        $options \\
        ${user_columns}

    if [[ ${out_option} != "--samout" ]]
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
