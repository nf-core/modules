process SLAMDUNK_ALL {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/slamdunk:0.4.3--py_0':
        'biocontainers/slamdunk:0.4.3--py_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)
    tuple val(meta4), path(filter_bed)

    output:
    tuple val(meta), path("outputs/map/*.bam")            , emit: bam
    tuple val(meta), path("outputs/filter/*.bam")         , emit: filtered_bam
    tuple val(meta), path("outputs/filter/*.bam.bai")     , emit: filtered_bai
    tuple val(meta), path("outputs/snp/*.vcf")            , emit: vcf, optional: true
    tuple val(meta), path("outputs/count/*.tsv")          , emit: tsv
    tuple val(meta), path("outputs/count/*_plus.bedgraph"), emit: plus_bedgraph
    tuple val(meta), path("outputs/count/*_mins.bedgraph"), emit: mins_bedgraph
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def filterbed = filter_bed ? "-fb ${filter_bed}" : ""
    """
    slamdunk \\
        all \\
        -r $fasta \\
        -b $bed \\
        -t $task.cpus \\
        -o outputs \\
        $args \\
        $filterbed \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slamdunk: \$(echo \$(slamdunk --version) | sed 's/^slamdunk //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p outputs/map
    mkdir -p outputs/filter
    mkdir -p outputs/snp
    mkdir -p outputs/count

    touch outputs/map/${prefix}.bam
    touch outputs/filter/${prefix}_filtered.bam
    touch outputs/filter/${prefix}_filtered.bam.bai
    touch outputs/snp/${prefix}.vcf
    touch outputs/count/${prefix}.tsv
    touch outputs/count/${prefix}_plus.bedgraph
    touch outputs/count/${prefix}_mins.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slamdunk: \$(echo \$(slamdunk --version) | sed 's/^slamdunk //')
    END_VERSIONS
    """
}
