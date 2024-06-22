process RIBOTISH_PREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotish:0.2.7--pyhdfd78af_0':
        'biocontainers/ribotish:0.2.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam_ribo), path(bai_ribo)
    tuple val(meta2), path(bam_ti), path(bai_ti)
    tuple val(meta3), path(fasta), path(gtf)
    tuple val(meta4), path(candidate_orfs)
    tuple val(meta5), path(para_ribo)
    tuple val(meta6), path(para_ti)

    output:
    tuple val(meta), path("*_pred.txt")        , emit: predictions
    tuple val(meta), path("*_all.txt")         , emit: all
    tuple val(meta), path("*_transprofile.py") , emit: transprofile
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    ribo_bam_cmd = ''
    ti_bam_cmd = ''
    if (bam_ribo){
        ribo_bam_cmd = "-b ${bam_ribo.join(',')}"
        if (para_ribo){
            ribo_bam_cmd += " --ribopara ${para_ribo.join(',')}"
        }
    }
    if (bam_ti){
        ti_bam_cmd = "-t ${bam_ti.join(',')}"
        if (para_tis){
            ti_bam_cmd += " --tisparapara  ${para_ti.join(',')}"
        }
    }
    """
    ribotish predict \\
        $ribo_bam_cmd \\
        $ti_bam_cmd \\
        -f $fasta \\
        -g $gtf \\
        -o ${prefix}_pred.txt \\
        --allresult ${prefix}_all.txt \\
        --transprofile ${prefix}_transprofile.py \\
        -p $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pred.txt
    touch ${prefix}_all.txt
    touch ${prefix}_transprofile.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | sed 's/ribotish //')
    END_VERSIONS
    """
}
