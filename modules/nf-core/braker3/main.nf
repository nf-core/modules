process BRAKER3 {
    tag "${meta.id}"
    label 'process_high'

    // Re. Conda from the BRAKER team:
    // Warning: installing GeneMark-ETP for BRAKER in conda environments has lead to multiple problems reported by users (Issues!).
    // We can not offer support for conda installations. Please use the singularity image instead.
    container "docker.io/teambraker/braker3:v3.0.7.5"

    input:
    tuple val(meta), path(fasta)
    path bam
    path rnaseq_sets_dirs
    path rnaseq_sets_ids
    path proteins
    path hintsfile

    output:
    tuple val(meta), path("$prefix/braker.gtf")         , emit: gtf
    tuple val(meta), path("$prefix/braker.codingseq")   , emit: cds
    tuple val(meta), path("$prefix/braker.aa")          , emit: aa
    tuple val(meta), path("$prefix/braker.log")         , emit: log
    tuple val(meta), path("$prefix/hintsfile.gff")      , emit: hintsfile   , optional: true
    tuple val(meta), path("$prefix/braker.gff3")        , emit: gff3        , optional: true
    tuple val(meta), path("$prefix/what-to-cite.txt")   , emit: citations
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args             ?: ''
    prefix          = task.ext.prefix           ?: "${meta.id}"
    def rna_ids     = rnaseq_sets_ids           ? "--rnaseq_sets_ids=$rnaseq_sets_ids"      : ''
    def rna_dirs    = rnaseq_sets_dirs          ? "--rnaseq_sets_dirs=$rnaseq_sets_dirs"    : ''
    def bam_arg     = bam                       ? "--bam=$bam"                              : ''
    def prot_arg    = proteins                  ? "--prot_seq=$proteins"                    : ''
    def hints       = hintsfile                 ? "--hints=$hintsfile"                      : ''
    def new_species = args.contains('--species')? ''                                        : '--species new_species'
    """
    cp -r \$AUGUSTUS_CONFIG_PATH \\
        augustus_config

    chmod -R a+w \\
        augustus_config

    perl -p -e 's/^(>\\S+).*\$/\$1/' \\
        $fasta \\
        > ${prefix}.name.only.genome.masked.fasta

    braker.pl \\
        --genome ${prefix}.name.only.genome.masked.fasta \\
        $new_species \\
        --workingdir $prefix \\
        --AUGUSTUS_CONFIG_PATH "\$(pwd)/augustus_config" \\
        --threads $task.cpus \\
        $rna_ids \\
        $rna_dirs \\
        $bam_arg \\
        $prot_arg \\
        $hints \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker3: \$(braker.pl --version 2>/dev/null | sed 's/braker.pl version //')
        augustus: \$(augustus --version |& sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p')
        genemark-etp: \$(echo "\$(gmetp.pl || echo '')" | sed -n 's/ETP version \\(.*\\)/\\1/p')
        prothint: \$(prothint.py --version | sed 's/prothint.py //1')
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args                         ?: ''
    prefix          = task.ext.prefix                       ?: "${meta.id}"
    def rna_ids     = rnaseq_sets_ids                       ? "--rnaseq_sets_ids=$rnaseq_sets_ids"      : ''
    def hints       = hintsfile                             ? "--hints=$hintsfile"                      : ''
    def touch_hints = (rna_ids || bam || proteins || hints) ? "touch $prefix/hintsfile.gff"             : ''
    def touch_gff   = args.contains('--gff3')               ? "touch $prefix/braker.gff3"               : ''
    """
    mkdir "$prefix"

    touch "$prefix/braker.gtf"
    touch "$prefix/braker.codingseq"
    touch "$prefix/braker.aa"
    $touch_hints
    touch "$prefix/braker.log"
    touch "$prefix/what-to-cite.txt"
    $touch_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker3: \$(braker.pl --version 2>/dev/null | sed 's/braker.pl version //')
        augustus: \$(augustus --version |& sed -n 's/AUGUSTUS (\\(.*\\)) is a gene .*/\\1/p')
        genemark-etp: \$(echo "\$(gmetp.pl || echo '')" | sed -n 's/ETP version \\(.*\\)/\\1/p')
        prothint: \$(prothint.py --version | sed 's/prothint.py //1')
    END_VERSIONS
    """
}
