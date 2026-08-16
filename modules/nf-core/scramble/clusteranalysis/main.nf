process SCRAMBLE_CLUSTERANALYSIS {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65d3a32dfd347b370e87589189717c75468e6d737b7cee6931e4dae21ce1a9cf/data':
        'community.wave.seqera.io/library/bioconductor-pwalign_scramble:31d27d3832b0689e' }"

    input:
    tuple val(meta) , path(clusters)
    tuple val(meta2), path(fasta)
    path mei_ref

    output:
    tuple val(meta), path("*_MEIs.txt")              , emit: meis_tab, optional:true
    tuple val(meta), path("*_PredictedDeletions.txt"), emit: dels_tab, optional:true
    tuple val(meta), path("*.vcf")                   , emit: vcf     , optional:true
    tuple val("${task.process}"), val('scramble'), val('1.0.2'), topic: versions, emit: versions_scramble
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def scramble_path = workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
        ? "\$CONDA_PREFIX/share/scramble" : "/opt/conda/share/scramble"

    def blastdb = args.contains("--eval-dels") ? "makeblastdb -in ${fasta} -parse_seqids -title ${fasta} -dbtype nucl -out ${fasta}" : ""
    def reference = fasta ? "--ref `pwd`/${fasta}" : ""

    // The default file for the MEI reference is a file that's inside the container
    def mei_reference = mei_ref ? "`pwd`${mei_ref}" : "${scramble_path}/resources/MEI_consensus_seqs.fa"

    """
    ${blastdb}

    Rscript --vanilla ${scramble_path}/bin/SCRAMble.R \\
        --install-dir ${scramble_path}/bin \\
        ${args} \\
        --cluster-file `pwd`/${clusters} \\
        ${reference} \\
        --mei-refs ${mei_reference} \\
        --out-name `pwd`/${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def create_deletions = args.contains("--eval-dels") ? "touch ${prefix}_PredictedDeletions.txt" : ""
    def create_vcf = fasta ? "touch ${prefix}.vcf ": ""
    """
    touch ${prefix}_MEIs.txt
    ${create_deletions}
    ${create_vcf}
    """
}
