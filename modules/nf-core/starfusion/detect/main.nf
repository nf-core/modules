process STARFUSION {
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
    tuple val(meta), path("*.fusion_predictions.tsv")                   , emit: fusions
    tuple val(meta), path("*.abridged.tsv")                             , emit: abridged
    tuple val(meta), path("*.coding_effect.tsv")     , optional: true   , emit: coding_effect
    path "versions.yml"                                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_arg = reads ? (meta.single_end ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}") : ""
    def junction_arg =  junction ? "-J ${junction}" : ""
    def args = task.ext.args ?: ''
    """
    STAR-Fusion \\
        --genome_lib_dir $reference \\
        $fastq_arg \\
        $junction_arg \\
        --CPU $task.cpus \\
        --examine_coding_effect \\
        --output_dir . \\
        $args

    mv star-fusion.fusion_predictions.tsv ${prefix}.starfusion.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${prefix}.starfusion.abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${prefix}.starfusion.abridged.coding_effect.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.starfusion.fusion_predictions.tsv
    touch ${prefix}.starfusion.abridged.tsv
    touch ${prefix}.starfusion.abridged.coding_effect.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}
