process TRAITAR_PFAMGET {
    tag "PFAM database"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/hmmer_prodigal_pandas_parallel_pruned:ccae2eabc2a54ac8' :
        'community.wave.seqera.io/library/hmmer_prodigal_pandas_parallel_pruned:ccae2eabc2a54ac8' }"

    output:
    path "pfam_data", emit: pfam_db
    tuple val("${task.process}"), val('traitar'), eval('traitar --version 2>&1 | tail -1'), topic: versions, emit: versions_traitar

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    traitar pfam pfam_data

    """

    stub:
    """
    mkdir -p pfam_data
    touch pfam_data/Pfam-A.hmm
    """
}
