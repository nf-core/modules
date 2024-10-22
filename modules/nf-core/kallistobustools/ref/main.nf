process KALLISTOBUSTOOLS_REF {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2' :
        'biocontainers/kb-python:0.28.2--pyhdfd78af_2' }"

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

    stub:
    if (workflow_mode == "standard") {
        """
        touch kb_ref_out.idx \\
        touch t2g.txt \\
        touch cdna.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    } else {
        """
        touch kb_ref_out.idx \\
        touch t2g.txt \\
        touch cdna.fa
        touch intron.fa \\
        touch cdna_t2c.txt \\
        touch intron_t2c.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
        END_VERSIONS
        """
    }
}
