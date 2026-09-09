process CHROMAP_CHROMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d39e0b3f00c5469ffc2ceef4bd76959fba6313064ed2408dd4ccac498022ad6/data' :
        'community.wave.seqera.io/library/chromap_samtools:b975c17adf0096ba' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    path barcodes
    path whitelist
    path chr_order
    path pairs_chr_order

    output:
    tuple val(meta), path("*.bed.gz")     , optional:true, emit: bed
    tuple val(meta), path("*.bam")        , optional:true, emit: bam
    tuple val(meta), path("*.tagAlign.gz"), optional:true, emit: tagAlign
    tuple val(meta), path("*.pairs.gz")   , optional:true, emit: pairs
    tuple val("${task.process}"), val('chromap'), eval("chromap --version 2>&1"), topic: versions, emit: versions_chromap
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()

    def file_extension = args.contains("--SAM") ? 'sam' : args.contains("--TagAlign")? 'tagAlign' : args.contains("--pairs")? 'pairs' : 'bed'
    if (barcodes) {
        args_list << "-b ${barcodes.join(',')}"
        if (whitelist) {
            args_list << "--barcode-whitelist $whitelist"
        }
    }
    if (chr_order) {
        args_list << "--chr-order $chr_order"
    }
    if (pairs_chr_order){
        args_list << "--pairs-natural-chr-order $pairs_chr_order"
    }
    def final_args = args_list.join(' ')
    def compression_cmds = "gzip -n ${prefix}.${file_extension}"
    if (args.contains("--SAM")) {
        compression_cmds = """
        samtools view $args2 -@ $task.cpus -bh \\
            -o ${prefix}.bam ${prefix}.${file_extension}
        rm ${prefix}.${file_extension}
        """
    }
    if (meta.single_end) {
        """
        chromap \\
            $final_args \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads.join(',')} \\
            -o ${prefix}.${file_extension}

        $compression_cmds
        """
    } else {
        """
        chromap \\
            $final_args \\
            -t $task.cpus \\
            -x $index \\
            -r $fasta \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${prefix}.${file_extension}

        $compression_cmds
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bed.gz
    touch ${prefix}.bam
    echo "" | gzip > ${prefix}.tagAlign.gz
    echo "" | gzip > ${prefix}.pairs.gz
    """
}
