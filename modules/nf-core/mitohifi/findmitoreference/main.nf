process MITOHIFI_FINDMITOREFERENCE {
    tag '$species'
    label 'process_low'

    // MitoHifi does not exist as conda package
    // A Docker image is available at the MitoHifi Github repository
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/marcelauliano/mitohifi:master':
        'ghcr.io/marcelauliano/mitohifi:master' }"

    input:
    val species
    val email
    val min_length

    output:
    path "*.fasta",                 emit: fasta
    path "*.gb",                    emit: gb
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    findMitoReference.py \\
        --species $species \\
        --email $email \\
        --min_length $min_length \\
        --outfolder .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.fasta
    touch ${prefix}.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
    """
}
