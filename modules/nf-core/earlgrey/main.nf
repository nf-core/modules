process EARLGREY {
    tag "${meta.id}"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/earlgrey:7.0.3--hd63eeec_0'
        : 'biocontainers/earlgrey:7.0.3--hd63eeec_0'}"

    input:
    tuple val(meta), path(genome)
    val famdb

    output:
    tuple val(meta), path("${meta.id}"), emit: output_directory
    //tuple val(meta), path("${meta.id}/*.fasta"), emit: fasta_files
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*.gff"), emit: gff
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*family*.kable"), emit: summary_family_kable
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*family*.txt"), emit: summary_family_txt
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*highLevel*.kable"), emit: summary_highLevel_kable
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*highLevel*.txt"), emit: summary_highLevel_txt
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*.summaryPie.pdf"), emit: summary_plot
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*classification*.pdf"), emit: classification_landscape_plot
    tuple val(meta), path("*/*/${meta.id}_summaryFiles/*split_class*.pdf"), emit: split_class_landscape_plot
    tuple val("${task.process}"), val('earlgrey'), eval('earlGrey -h | grep version | cut -d" " -f3'), topic: versions, emit: versions_earlgrey

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = genome.toString() - ~/\.gz$/
    def decompress = genome.getExtension() == "gz" ? "gunzip -c ${genome} > ${input}" : ""
    def famdb_dir = "/usr/local/share/RepeatMasker/Libraries/famdb/"
    def localdb = famdb ? "mv */*.h5.gz ${famdb_dir}" : ""

    """
    # unzip genome
    ${decompress}

    # Move partitions ot famdb directory
    ${localdb}

    # capture current working directory
    INPUT_DIR=\$PWD

    # Change directory to the famdb library location
    cd ${localdb}

    # download the required partitions from Dfam 3.9
    # curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0].h5.gz'

    # gunzip -f *.gz

    # move up to RepeatMasker main directory
    cd /usr/local/share/RepeatMasker/

    # Rerun RepeatMasker configuration
    perl ./configure \\
            -libdir /usr/local/share/RepeatMasker/Libraries/ \\
            -trf_prgm /usr/local/bin/trf \\
            -rmblast_dir /usr/local/bin \\
            -hmmer_dir /usr/local/bin \\
            -abblast_dir /usr/local/bin \\
            -crossmatch_dir /usr/local/bin \\
            -default_search_engine rmblast

    cd \$INPUT_DIR

    mkdir ${prefix}

    set +e
    earlGrey \\
        ${args} \\
        -g ${input} \\
        -s ${prefix} \\
        -t ${task.cpus} \\
        -o ${prefix} > earlgrey.log 2>&1
    EXIT_STATUS=\$?
    set -e

    # Check for specific message in the log and run again (necessary for 0 partition)
    if [ \$EXIT_STATUS -eq 2 ] && grep -q "If you only want partition 0" earlgrey.log; then
        echo "Message found in first run, executing second run..."
        earlGrey \\
            ${args} \\
            -g ${input} \\
            -s ${prefix} \\
            -t ${task.cpus} \\
            -o ${prefix}
    else
        echo "Message not found, skipping second run"
    fi

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo ${args}

    mkdir ${prefix}
    """
}
