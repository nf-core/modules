process PARABRICKS_STAR {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(index)
    tuple val(meta3), path(genome_lib_dir)

    output:
    tuple val(meta), path("*.bam")                  , emit: bam              , optional:true
    tuple val(meta), path("*.bai")                  , emit: bai              , optional:true
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"

    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''

    // Mostly for testing purposes, to allow the genome directory to be passed as tarball
    if genome_lib_dir.endswith('tar.gz') || genome_lib_dir.endswith('tar') || genome_lib_dir.endswith('tgz') || genome_lib_dir.endswith('zip') {
        genome_lib_dir = "genome_lib_dir"
        """
        mkdir -p $genome_lib_dir
        tar -xzf ${genome_lib_dir} -C $genome_lib_dir
        """
    } 

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp $fasta \$INDEX

    pbrun \\
        rna_fq2bam  \\
        --ref \$INDEX \\
        $in_fq_command \\
        --output-dir $prefix \\
        --genome-lib-dir ${genome_lib_dir} \\
        --out-bam ${prefix}.bam \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    $qc_metrics_output
    $duplicate_metrics_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
