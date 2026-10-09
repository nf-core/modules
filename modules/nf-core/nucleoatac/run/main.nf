process NUCLEOATAC_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nucleoatac:1.0.0--py310h3479294_0':
        'quay.io/biocontainers/nucleoatac:1.0.0--py310h3479294_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bed), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.nucpos.bed*")                       , emit: nucpos, optional: true
    tuple val(meta), path("*.nucpos.redundant.bed*")             , emit: nucpos_redundant, optional: true
    tuple val(meta), path("*.nfrpos.bed*")                       , emit: nfrpos, optional: true
    tuple val(meta), path("*.nucmap_combined.bed*")              , emit: nucmap_combined, optional: true
    tuple val(meta), path("*.occpeaks.bed*")                     , emit: occpeaks, optional: true
    tuple val(meta), path("*.occ.bedgraph*")                     , emit: occupancy, optional: true
    tuple val(meta), path("*.occ.lower_bound.bedgraph*")         , emit: occupancy_lower_bound, optional: true
    tuple val(meta), path("*.occ.upper_bound.bedgraph*")         , emit: occupancy_upper_bound, optional: true
    tuple val(meta), path("*.nucleoatac_signal.bedgraph*")       , emit: signal, optional: true
    tuple val(meta), path("*.nucleoatac_signal.smooth.bedgraph*"), emit: signal_smooth, optional: true
    tuple val(meta), path("*.fragmentsizes.txt")                 , emit: fragment_sizes, optional: true
    tuple val(meta), path("*.nuc_dist.txt")                      , emit: nuc_dist, optional: true
    tuple val(meta), path("*.VMat")                              , emit: vmat, optional: true
    tuple val(meta), path("*.eps")                               , emit: plots, optional: true
    tuple val("${task.process}"), val('nucleoatac'), eval('nucleoatac --version | sed "s/^nucleoatac //"'), topic: versions, emit: versions_nucleoatac

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_line = args ? " \\\n        ${args}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ln -s $bam ${prefix}.bam
    if [[ "$bai" == *.csi ]]; then
        ln -s $bai ${prefix}.bam.csi
    else
        ln -s $bai ${prefix}.bam.bai
    fi

    nucleoatac run \\
        --bed $bed \\
        --bam ${prefix}.bam \\
        --fasta $fasta \\
        --out ${prefix} \\
        --cores ${task.cpus}${args_line}

    for file in ${prefix}*.bed ${prefix}*.bedgraph; do
        if [[ -f "\$file" ]]; then
            gzip -f "\$file"
        fi
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.nucpos.bed.gz
    echo "" | gzip > ${prefix}.nucpos.redundant.bed.gz
    echo "" | gzip > ${prefix}.nfrpos.bed.gz
    echo "" | gzip > ${prefix}.nucmap_combined.bed.gz
    echo "" | gzip > ${prefix}.occpeaks.bed.gz
    echo "" | gzip > ${prefix}.occ.bedgraph.gz
    echo "" | gzip > ${prefix}.occ.lower_bound.bedgraph.gz
    echo "" | gzip > ${prefix}.occ.upper_bound.bedgraph.gz
    echo "" | gzip > ${prefix}.nucleoatac_signal.bedgraph.gz
    echo "" | gzip > ${prefix}.nucleoatac_signal.smooth.bedgraph.gz
    printf "insert_size\\tcount\\n150\\t1\\n" > ${prefix}.fragmentsizes.txt
    printf "insert_size\\tdensity\\n150\\t1\\n" > ${prefix}.nuc_dist.txt
    printf "position\\tsignal\\n0\\t0\\n" > ${prefix}.VMat
    touch ${prefix}.nuc_dist.eps ${prefix}.occ_fit.eps
    """
}
