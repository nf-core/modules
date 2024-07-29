process SENTIEON_READWRITER {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sentieon:202308.02--ffce1b7074ce9924' :
        'nf-core/sentieon:202308.02--c641bc397cbf79d5' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.${format}"),                     emit: output
    tuple val(meta), path("*.${index}") ,                     emit: index
    tuple val(meta), path("*.${format}"), path("*.${index}"), emit: output_index
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:


    def args            = task.ext.args   ?: ''
    def args2           = task.ext.args2  ?: ''
    def input_str       = input.sort().collect{"-i $it"}.join(' ')
    def reference       = fasta ? "-r $fasta" : ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    format              = input.extension == "bam" ? "bam"     : "cram"
    index               = format          == "bam" ? "bam.bai" : "cram.crai"
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
        ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:


    def prefix = task.ext.prefix ?: "${meta.id}"
    format = input.extension == "bam" ? "bam"     : "cram"
    index  = format          == "bam" ? "bam.bai" : "cram.crai"
    """

    touch ${prefix}.${format}
    touch ${prefix}.${index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
