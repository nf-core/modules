/**
This module describes the custom Rastair process for
parsing the rastair call output
and converting into MethylKit and Bismark digestible formats.
*/

process CONVERT_TO_METHYLKIT {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f0a47f3c0c4f521ed0623cd709c51fb1ece4df1fb4bd85c75d04e0383a8c5d4/data' :
        'community.wave.seqera.io/library/rastair:0.8.2--b09a8e25a0d53059' }"

    input:
    tuple val(meta), path(rastair_call_txt)

    output:
    tuple val(meta), path("*methylkit.txt.gz"), emit: methylkit
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat ${rastair_call_txt} | /opt/conda/share/rastair/scripts/rastair_call_to_methylkit.sh | gzip -c > ${prefix}.rastair_methylkit.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
