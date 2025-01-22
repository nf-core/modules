process PARABRICKS_DEEPVARIANT {
    tag "$meta.id"
    label 'gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1"

    input:
    tuple val(meta), path(bam), path(bai), path(target_region_bed)
    tuple path(fasta), path(fai)
    path model_file
    path proposed_variants


    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix     = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def output_file = args =~ "gvcf" ? "${prefix}.genome.vcf" : "${prefix}.vcf"
    def target_region_bed_command = target_region_bed ? target_region_bed.collect{"--interval-file $it"}.join(' ') : ""
    def proposed_variants_option = proposed_variants ? "--proposed-variants $proposed_variants" : ""
    def model_file_option = model_file ? "--pb-model-file $model_file" : ""
    def wes_model_option = (params.assay == "WES") ? "--use-wes-model" : ""

    """

    logfile=run.log
    exec > >(tee \$logfile)
    exec 2>&1

    echo "pbrun deepvariant --ref $fasta --in-bam $bam --out-variants $output_file --num-gpus $task.accelerator.request $wes_model_option $target_region_bed_command $proposed_variants_option $model_file_option $args"

    pbrun \\
        deepvariant \\
        --ref $fasta \\
        --in-bam $bam \\
        --out-variants $output_file \\
        --num-gpus $task.accelerator.request \\
        $wes_model_option \\
        $target_region_bed_command \\
        $proposed_variants_option \\
        $model_file_option \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = args =~ "gvcf" ? "${prefix}.genome.vcf" : "${prefix}.vcf"
    """
    touch run.log
    touch $output_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}