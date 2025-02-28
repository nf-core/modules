process SENTIEON_READWRITER {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a6/a64461f38d76bebea8e21441079e76e663e1168b0c59dafee6ee58440ad8c8ac/data' :
        'community.wave.seqera.io/library/sentieon:202308.03--59589f002351c221' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}"),                             emit: output
    tuple val(meta), path("${prefix}.${index}"),                    emit: index
    tuple val(meta), path("${prefix}"), path("${prefix}.${index}"), emit: output_index
    path  "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:


    def args            = task.ext.args   ?: ''
    def args2           = task.ext.args2  ?: ''
    def input_str       = input.sort{ it.getName() }.collect{"-i $it"}.join(' ')
    def reference       = fasta ? "-r $fasta" : ''

    // bam -> bam: prefix = "<filename>.bam"
    // bam -> cram: prefix = "<filename>.cram"
    // cram -> cram: prefix = "<filename>.cram"
    prefix          = task.ext.prefix ?: "${meta.id}.bam"
    index           = prefix.tokenize('.')[-1] == "bam" ? "bai" : "crai"

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""
    """
    $sentieonLicense

    sentieon \\
        driver \\
        -t $task.cpus \\
        $reference \\
        $args \\
        $input_str \\
        --algo ReadWriter \\
        $args2 \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix ?: "${meta.id}.cram"
    index           = prefix.tokenize('.')[-1] == "bam" ? "bai" : "crai"
    """

    touch ${prefix}
    touch ${prefix}.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
