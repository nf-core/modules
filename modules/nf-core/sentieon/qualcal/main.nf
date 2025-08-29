process SENTIEON_QUALCAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/20/2050c5321a5426e31b9ed1e3e98356913fe3c316a7ef02c4fb872983a730db6f/data'
        : 'community.wave.seqera.io/library/sentieon_gnuplot:a6da525a6c9ce6e3'}"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(known_sites)
    tuple val(meta5), path(known_sites_tbi)
    tuple val(meta6), path(recalibration_table)
    val generate_recalibrated_bams

    output:
    tuple val(meta), path("*.table"),      emit: table, optional: true
    tuple val(meta), path("*.table.post"), emit: table_post, optional: true
    tuple val(meta), path("*.{cram,bam}"), emit: recal_alignment, optional: true
    tuple val(meta), path("*.csv"),        emit: csv, optional: true
    tuple val(meta), path("*.pdf"),        emit: pdf, optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = input.collect { "-i ${it}" }.join(' ')
    def knownSites = known_sites ? known_sites.collect { "-k ${it}" }.join(' ') : ""
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""

    // Create recalibration table. Actual base quality recalibration can be done during Variant calling with Sentieon
    if (!recalibration_table) {
        """
        ${sentieonLicense}

        sentieon driver \\
            -t ${task.cpus} \\
            -r ${fasta} \\
            ${input_list} \\
            --algo QualCal \\
            ${args} \\
            ${knownSites} \\
            ${prefix}.table

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
    }
    else {
        // Runs basequality recalibration with a previous generated table
        def file_suffix = input.name.endsWith(".bam") ? "bam" : "cram"
        def recalibrated_bam = generate_recalibrated_bams ? "--algo ReadWriter ${prefix}.recalibrated.${file_suffix}" : ""
        """
        ${sentieonLicense}

        sentieon driver \\
            -t ${task.cpus} \\
            -r ${fasta} \\
            -i ${input} \\
            -q ${recalibration_table} \\
            --algo QualCal \\
            ${args} \\
            ${knownSites} \\
            ${prefix}.table.post \\
            ${recalibrated_bam} \\

        sentieon driver \\
            -t ${task.cpus} \\
            --algo QualCal \\
            --plot \\
            --before ${prefix}.table \\
            --after ${prefix}.table.post \\
            ${args} \\
            ${prefix}.csv

        sentieon plot QualCal -o ${prefix}.pdf ${prefix}.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_suffix = input.name.endsWith(".bam") ? "bam" : "cram"
    def recalibrated_bam = generate_recalibrated_bams ? "${prefix}.recalibrated.${file_suffix}" : ""
    def recalibrated_bam_cmd = generate_recalibrated_bams ? "touch ${recalibrated_bam}" : ""
    """
    touch ${prefix}.table
    touch ${prefix}.table.post
    ${recalibrated_bam_cmd}
    touch ${prefix}.csv
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
