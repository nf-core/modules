process STARFUSION_DETECT {
    tag "$meta.id"
    label 'process_high'
// WARN: Tool is reporting wrong version number in current release, please check if this is fixed when updating.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75d085bf2a8e40c6693b357800eef0f9568f661226d0888339bc77f7852234bb/data' :
        'community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:e285bb3eb373b9a7'}"

    input:
    tuple val(meta), path(reads), path(junction)
    path reference

    output:
    tuple val(meta), path("*.fusion_predictions.tsv"), emit: fusions
    tuple val(meta), path("*.abridged.tsv")          , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , emit: coding_effect, optional: true
    // STAR-Fusion's own `--version` prints '1.15.0' even for the 1.15.1 release, so
    // we emit the correct pinned version (from environment.yml) rather than the
    // tool's self-reported string.
    tuple val("${task.process}"), val('star-fusion'), eval("echo 1.15.1"), topic: versions, emit: versions_starfusion

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}.starfusion"
    def fastq_arg    = reads ? (meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}") : ""
    def junction_arg =  junction ? "-J ${junction}" : ""
    def args         = task.ext.args ?: ''
    """
    STAR-Fusion \\
        --genome_lib_dir $reference \\
        $fastq_arg \\
        $junction_arg \\
        --CPU $task.cpus \\
        --output_dir . \\
        $args

    mv star-fusion.fusion_predictions.tsv ${prefix}.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${prefix}.abridged.tsv
    if [ -f star-fusion.fusion_predictions.abridged.coding_effect.tsv ]; then
        mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${prefix}.abridged.coding_effect.tsv
    fi
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}.starfusion"
    """
    touch ${prefix}.fusion_predictions.tsv
    touch ${prefix}.abridged.tsv
    touch ${prefix}.abridged.coding_effect.tsv
    """
}
