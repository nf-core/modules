process PARAPHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0':
        'biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0' }"

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
    touch ${prefix}_paraphase_vcfs/${prefix}_stub.vcf.gz.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        paraphase: \$(paraphase --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
