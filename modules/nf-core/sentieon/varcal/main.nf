process SENTIEON_VARCAL {
    tag "${meta.id}"
    label 'process_low'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/80/80ccb05eb4f1a193a3bd99c4da90f55f74ea6556c25f154e53e1ff5a6caa372d/data'
        : 'community.wave.seqera.io/library/sentieon:202503--5e378058d837c58c'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path resource_vcf
    path resource_tbi
    val labels
    path fasta
    path fai

    output:
    tuple val(meta), path("*.recal"), emit: recal
    tuple val(meta), path("*.idx"), emit: idx
    tuple val(meta), path("*.tranches"), emit: tranches
    tuple val(meta), path("*plots.R"), emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Process labels to create the command string
    // labels is a list. Here is an example of what labels might look like:
    // ['--resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz', '--resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz']
    def labels_command = ''
    def labels_input = labels
    if (labels_input instanceof String && !labels_input.trim().isEmpty()) {
        // Process string input
        def resourceStrings = labels_input.split('--resource:').findAll()
        def processedResources = resourceStrings.collect { resource_string ->
            def items = resource_string.split(' ', 2)
            if (items.size() != 2) {
                error("Expected the resource string '${resource_string}' to contain two elements separated by a space.")
            }
            "--resource ${items[1]} --resource_param ${items[0].replaceFirst('^--resource:', '')}"
        }
        labels_command = processedResources.join(' ')
    }
    else if (labels_input instanceof List) {
        // Process list input
        def processedResources = labels_input.collect { label ->
            def cleanedLabel = label.replaceFirst('^--resource:', '')
            def items = cleanedLabel.split(' ', 2)
            if (items.size() != 2) {
                error("Expected the resource string '${cleanedLabel}' to contain two elements separated by a space.")
            }
            "--resource ${items[1]} --resource_param ${items[0].replaceFirst('^--resource:', '')}"
        }
        labels_command = processedResources.join(' ')
    }
    else if (labels_input != null) {
        error("Expected 'labels' to be either a String or a List, but got ${labels_input.getClass()}")
    }

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver -r ${fasta}  --algo VarCal \\
        -v ${vcf} \\
        --tranches_file ${prefix}.tranches \\
        ${labels_command} \\
        ${args} \\
        ${prefix}.recal

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal
    touch ${prefix}.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
