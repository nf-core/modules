process PARAPHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9f4f65b0d84c9cc1a7f62830a13170cf46293f264ee27ded601346a01befa58/data':
        'community.wave.seqera.io/library/minimap2_paraphase_samtools:be7ebf7454267f19' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(config)

    output:
    tuple val(meta), path("*.paraphase.json")                           , emit: json
    tuple val(meta), path("*.paraphase.bam")                            , emit: bam
    tuple val(meta), path("*.paraphase.bam.bai")                        , emit: bai
    tuple val(meta), path("${prefix}_paraphase_vcfs/*.vcf.gz")          , emit: vcf      , optional: true
    tuple val(meta), path("${prefix}_paraphase_vcfs/*.vcf.gz.{csi,tbi}"), emit: vcf_index, optional: true
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def config_file = config ? "--config $config" : ""
    """
    paraphase \\
        $args \\
        --threads $task.cpus \\
        --bam $bam \\
        --reference $fasta \\
        --prefix $prefix \\
        $config_file \\
        --out .

    for vcf in ${prefix}_paraphase_vcfs/*.vcf;
    do
        bgzip \\
            $args2 \\
            --threads $task.cpus \\
            \$vcf;
        tabix \\
            $args3 \\
            --threads $task.cpus \\
            \$vcf.gz;
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        paraphase: \$(paraphase --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args3 = task.ext.args3 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def index = args3.contains('--csi') ? 'csi' : 'tbi'
    """
    mkdir ${prefix}_paraphase_vcfs

    touch ${prefix}.paraphase.json
    touch ${prefix}.paraphase.bam
    touch ${prefix}.paraphase.bam.bai
    echo '' | gzip > ${prefix}_paraphase_vcfs/${prefix}_stub.vcf.gz
    echo "" | gzip > ${prefix}_paraphase_vcfs/${prefix}_stub.vcf.gz.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        paraphase: \$(paraphase --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
