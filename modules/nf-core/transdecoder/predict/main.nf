process TRANSDECODER_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.7.1--pl5321hdfd78af_0':
        'quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(fold)

    output:
    tuple val(meta), path("*.transdecoder.pep")  , emit: pep
    tuple val(meta), path("*.transdecoder.gff3") , emit: gff3
    tuple val(meta), path("*.transdecoder.cds")  , emit: cds
    tuple val(meta), path("*.transdecoder.bed")  , emit: bed
    tuple val("${task.process}"), val('transdecoder'), eval("TransDecoder.Predict --version | sed 's/TransDecoder.Predict //'"), emit: versions_transdecoder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # \$fold is staged into TRANSDECODER_LONGORF's own output -- a symlink under local/shared-filesystem
    # staging, but a real copy under stageInMode 'copy' (the default for cloud storage without Fusion).
    # TransDecoder.Predict writes its checkpoints and intermediates there, which mutates that other
    # task's output and makes this task uncacheable across -resume (nf-core/modules#12799). Give it a
    # private, writable directory of symlinks instead. mv handles both staging modes: it moves a
    # symlink as a symlink, and renames a real directory in place, so nothing is deleted either way.
    shopt -s dotglob nullglob
    mv "$fold" "${fold}.staged"
    real_fold=\$(readlink -f "${fold}.staged")
    [ -n "\$(ls -A "\$real_fold")" ] || { echo "TRANSDECODER_LONGORF produced an empty directory: \$real_fold" >&2; exit 1; }
    mkdir "$fold"
    ln -s "\$real_fold"/* "$fold"/
    shopt -u dotglob nullglob

    TransDecoder.Predict \\
        $args \\
        -O . \\
        -t \\
        $fasta
    """

    stub:
    def fasta_no_gz = fasta.toString() - '.gz'
    """
    touch ${fasta_no_gz}.transdecoder.pep
    touch ${fasta_no_gz}.transdecoder.gff3
    touch ${fasta_no_gz}.transdecoder.cds
    touch ${fasta_no_gz}.transdecoder.bed
    """
}
