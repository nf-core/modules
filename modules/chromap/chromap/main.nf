process CHROMAP_CHROMAP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::chromap=0.2.1 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:bd74d08a359024829a7aec1638a28607bbcd8a58-0' :
        'quay.io/biocontainers/mulled-v2-1f09f39f20b1c4ee36581dc81cc323c70e661633:bd74d08a359024829a7aec1638a28607bbcd8a58-0' }"

    input:
    tuple val(meta), path(reads)
    path fasta
    path index
    path barcodes
    path whitelist
    path chr_order
    path pairs_chr_order

    output:
    tuple val(meta), path("*.bed.gz")     , optional:true, emit: bed
    tuple val(meta), path("*.bam")        , optional:true, emit: bam
    tuple val(meta), path("*.tagAlign.gz"), optional:true, emit: tagAlign
    tuple val(meta), path("*.pairs.gz")   , optional:true, emit: pairs
    path "versions.yml"                                  , emit: versions

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            chromap: \$(echo \$(chromap --version 2>&1))
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            chromap: \$(echo \$(chromap --version 2>&1))
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
