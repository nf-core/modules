process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sentieon:202308.02--ffce1b7074ce9924' :
        'nf-core/sentieon:202308.02--c641bc397cbf79d5' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam_and_bai
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    if [ -n "\${SENTIEON_LICENSE_BASE64}" ]; then
        # NOTE: This is how pipeline users will test out Sentieon
        echo "Localhost license file"
        # The license file is stored as a nextflow variable like, for instance, this:
        # nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
        export SENTIEON_LICENSE=\$(mktemp)
        echo -e "\$SENTIEON_LICENSE_BASE64" | base64 -d > \$SENTIEON_LICENSE
    fi


    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    sentieon bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | sentieon util sort -r $fasta -t $task.cpus -o ${prefix}.bam --sam2bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
