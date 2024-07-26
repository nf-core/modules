process SENTIEON_QUALCAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/sentieon:202308.02--ffce1b7074ce9924' :
        'nf-core/sentieon:202308.02--c641bc397cbf79d5' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(known_sites) //optional
    tuple val(meta4), path(recal_table) // not optional?!
    tuple val(meta5), path(recal_table_post) // not optional?!
    val generate_recalibrated_bams // false?

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true  //recalibrated bam: output is optional
    tuple val(meta), path("*.report"), emit: report, optional: true
    tuple val(meta), path ("recal_table"), emit: recal_table, optional: true 
    tuple val(meta), path ("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known_sites = known_sites ? known_sites.collect{"-k $it"}.join(' ') : ""
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64 ?
        "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; " :
        ""

    // Temprorary output files:
    // recal_table_post
    // recal_csv

    // Actual base quality recalibration can be done during Variant calling with Sentieon
    if() {
        """
        $sentieonLicense

        sentieon driver --algo QualCal \\
            $args \\
            -t $task.cpus \\
            -r $fasta \\
            -i $bam \\
            $known_sites \\
            ${prefix}.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
    } else {
    // Runs basequality recalibration in one step
    //TODO bam and cram should both work, this is also still optional
        def recalibrated_bam = generate_reclibrated_bams ? "--algo ReadWriter ${prefix}.recalibrated.cram" : ""
        """
        $sentieonLicense

        sentieon driver --algo QualCal \\
            $args \\
            -t $task.cpus \\
            -r $fasta \\
            -i $bam \\
            $known_sites \\
            -q ${prefix}.table \\
            $recalibrated_bam \\
            ${prefix}.table.post

        sentieon driver --algo QualCal \\
            $args \\
            -t $task.cpus \\
            --plot \\
            --before ${prefix}.table \\
            --after ${prefix}.table.post
            ${prefix}.csv

        sentieon plot QualCal -o ${prefix}.pdf ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
