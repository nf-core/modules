process PARABRICKS_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy' // needed by the module to work properly can be removed when fixed upstream - Issue #7226

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(interval_file)
    tuple val(ref_meta), path(fasta)

    output:
    tuple val(meta), path("*.vcf"),      optional: true, emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), optional: true, emit: gvcf
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = ("gvcf" =~ task.ext.args)? "${prefix}.g.vcf.gz" : "${prefix}.vcf"
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    pbrun \\
        haplotypecaller \\
        --ref $fasta \\
        --in-bam $input \\
        --out-variants $output_file \\
        $interval_file_command \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = ("--gvcf" =~ task.ext.args)? "echo '' | gzip > ${prefix}.g.vcf.gz" : "touch ${prefix}.vcf"
    """
    $output_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
