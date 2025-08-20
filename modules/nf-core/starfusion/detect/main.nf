process STARFUSION_DETECT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75d085bf2a8e40c6693b357800eef0f9568f661226d0888339bc77f7852234bb/data' :
        'community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:e285bb3eb373b9a7'}"

    input:
    tuple val(meta), path(reads), path(junction)
    path reference

    output:
    tuple val(meta), path("*.fusion_predictions.tsv"), emit: fusions
    tuple val(meta), path("*.abridged.tsv")          , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , emit: coding_effect, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix       = task.ext.prefix ?: "${meta.id}.starfusion"
    def fastq_arg    = reads ? (meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}") : ""
    def junction_arg =  junction ? "-J ${junction}" : ""
    def args         = task.ext.args ?: ''
    def VERSION      = '1.15.1' // WARN: This is the actual version of the STAR-FUSION, but version information of tool is not updated and prints '1.15.0'
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}.starfusion"
    def VERSION = '1.15.1'
    """
    touch ${prefix}.fusion_predictions.tsv
    touch ${prefix}.abridged.tsv
    touch ${prefix}.abridged.coding_effect.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: $VERSION
    END_VERSIONS
    """
}
