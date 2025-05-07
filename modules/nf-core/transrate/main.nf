process TRANSRATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c7d34ed4d4e8b30bc4216c753fcf3084d3bab43f8f2e80eb63b6184fbd25f04/data':
        'community.wave.seqera.io/library/snap-aligner_transrate:501869af4f81472b' }"

    input:
    tuple val(meta) , path(fasta)     // assembly file
    tuple val(meta2), path(reads)     // processed reads
    tuple val(meta3), path(reference) // reference proteins or transcripts fasta

    output:
    tuple val(meta), path("${prefix}"), emit: transrate_results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def transrate_assembly = fasta ? "--assembly $fasta" : ''
    def transrate_reads = reads ? ( meta2.single_end ? "--left $reads" : "--left ${reads[0]} --right ${reads[1]}" ) : "" // TODO test with meta.single_end
    def transrate_reference = reference ? "--reference $reference" : ''
    """
    transrate \\
        $args \\
        $transrate_assembly \\
        $transrate_reads \\
        $transrate_reference \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: \$(transrate --version 2>/dev/null | tail -n1)
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //;s/.\$//')
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    touch ${prefix}/test.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: \$(transrate --version 2>/dev/null | tail -n1)
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //;s/.\$//')
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        blast: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """
}
