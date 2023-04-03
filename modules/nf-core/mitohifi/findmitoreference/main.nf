process MITOHIFI_FINDMITOREFERENCE {
    tag '$species'
    label 'process_low'

    // Docker image available at the biocontainers Dockerhub
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/mitohifi:3.0.0_cv1':
        'docker.io/biocontainers/mitohifi:3.0.0_cv1' }"

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
    touch accession.fasta
    touch accession.gb

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
    """
}
