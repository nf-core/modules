process RSEM_CALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::rsem=1.3.3 bioconda::star=2.7.6a" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0' :
        'quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    tuple val(meta), path("*.log")             , emit: logs
    path  "versions.yml"                       , emit: versions

    tuple val(meta), path("*.STAR.genome.bam")       , optional:true, emit: bam_star
    tuple val(meta), path("${prefix}.genome.bam")    , optional:true, emit: bam_genome
    tuple val(meta), path("${prefix}.transcript.bam"), optional:true, emit: bam_transcript

    script:
    def args = task.ext.args ?: ''
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--strandedness forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--strandedness reverse'
    }
    def paired_end = meta.single_end ? "" : "--paired-end"
    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/.grp//'`
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --temporary-folder ./tmp/ \\
        $strandedness \\
        $paired_end \\
        $args \\
        $reads \\
        \$INDEX \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
