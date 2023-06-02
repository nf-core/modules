process SCRAMBLE_CLUSTERANALYSIS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::scramble=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scramble:1.0.1--h779adbc_1':
        'biocontainers/scramble:1.0.1--h779adbc_1' }"

    input:
    tuple val(meta), path(clusters)
    path fasta
    path mei_ref

    output:
    tuple val(meta), path("*_MEIs.txt")                 , optional:true, emit: meis_tab
    tuple val(meta), path("*_PredictedDeletions.txt")   , optional:true, emit: dels_tab
    tuple val(meta), path("*.vcf")                      , optional:true, emit: vcf
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def blastdb = args.contains("--eval-dels") ? "makeblastdb -in ${fasta} -parse_seqids -title ${fasta} -dbtype nucl -out ${fasta}" : ""
    def reference = fasta ? "--ref `pwd`/${fasta}" : ""

    // The default file for the MEI reference is a file that's inside the container
    def mei_reference = mei_ref ? "`pwd`/${mei_ref}" : "/usr/local/share/scramble/resources/MEI_consensus_seqs.fa"

    def blastdb_version = args.contains("--eval-dels") ? "makeblastdb: \$(echo \$(makeblastdb -version 2>&1) | head -n 1 | sed 's/^makeblastdb: //; s/+ Package.*\$//')" : ""
    """
    ${blastdb}

    Rscript --vanilla /usr/local/share/scramble/bin/SCRAMble.R \\
        --install-dir /usr/local/share/scramble/bin \\
        ${args} \\
        --cluster-file `pwd`/${clusters} \\
        ${reference} \\
        --mei-refs ${mei_reference} \\
        --out-name `pwd`/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scramble: ${VERSION}
        ${blastdb_version}
    END_VERSIONS
    """
}
