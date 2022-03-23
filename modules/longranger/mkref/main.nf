
process LONGRANGER_MKREF {
    tag "$meta.id"
    label 'process_medium'

     if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using longranger"
    }
    if ( workflow.containerEngine == 'singularity' || \
            workflow.containerEngine == 'docker' ) {
        exit 1, "Longranger can not be run in container environment"
    }

    memory '100 GB'

    input:
    tuple val(meta), path(reference)

    output:
    tuple val(meta), path("refdata-*"),                          emit: folder
    tuple val(meta), path("refdata-*/fasta/genome.fa.sa"),       emit: sa
    tuple val(meta), path("refdata-*/fasta/genome.fa.amb"),      emit: amb
    tuple val(meta), path("refdata-*/fasta/genome.fa.ann"),      emit: ann
    tuple val(meta), path("refdata-*/fasta/genome.fa.bwt"),      emit: bwt
    tuple val(meta), path("refdata-*/fasta/genome.fa.fai"),      emit: fai
    tuple val(meta), path("refdata-*/fasta/genome.fa.flat"),     emit: flat
    tuple val(meta), path("refdata-*/fasta/genome.fa.gdx"),      emit: gdx
    tuple val(meta), path("refdata-*/fasta/genome.fa.pac"),      emit: pac
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longranger mkref $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    longranger align --version
    END_VERSIONS
    """
}
