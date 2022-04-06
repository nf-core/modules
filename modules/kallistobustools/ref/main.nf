process KALLISTOBUSTOOLS_REF {
    tag "$fasta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::kb-python=0.26.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.26.3--pyhdfd78af_0' :
        'quay.io/biocontainers/kb-python:0.26.3--pyhdfd78af_0' }"

    input:
    path fasta
    path gtf
    val  workflow_mode

    output:
    path "versions.yml"   , emit: versions
    path "kb_ref_out.idx" , emit: index
    path "t2g.txt"        , emit: t2g
    path "cdna.fa"        , emit: cdna
    path "intron.fa"      , optional:true, emit: intron
    path "cdna_t2c.txt"   , optional:true, emit: cdna_t2c
    path "intron_t2c.txt" , optional:true, emit: intron_t2c

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (workflow_mode == "standard") {
        """
        kb \\
            ref \\
            -i kb_ref_out.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            --workflow $workflow_mode \\
            $fasta \\
            $gtf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    } else {
        """
        kb \\
            ref \\
            -i kb_ref_out.idx \\
            -g t2g.txt \\
            -f1 cdna.fa \\
            -f2 intron.fa \\
            -c1 cdna_t2c.txt \\
            -c2 intron_t2c.txt \\
            --workflow $workflow_mode \\
            $fasta \\
            $gtf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    }
}
