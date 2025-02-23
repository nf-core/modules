process PARABRICKS_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy' // needed by the module to work properly can be removed when fixed upstream - Issue #7226

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(bam), path(bai), path(interval_file)
    tuple val(ref_meta), path(fasta)
    path model_file

    output:
    tuple val(meta), path("*.vcf"), optional: true, emit: vcf
    tuple val(meta), path("*.g.vcf"), optional: true, emit: gvcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix     = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def output_file = ("--gvcf" =~ task.ext.args)? "${prefix}.g.vcf" : "${prefix}.vcf"
    def interval_file_option = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def model_command = model_file ? "--pb-model-file $model_file" : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''

    """
    pbrun \\
        deepvariant \\
        --ref $fasta \\
        --in-bam $bam \\
        --out-variants $output_file \\
        ${interval_file_option} \\
        ${num_gpus} \\
        ${model_command} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = ("--gvcf" =~ task.ext.args)? "echo '' | gzip > ${prefix}.g.vcf" : "touch ${prefix}.vcf"
    """
    $output_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
