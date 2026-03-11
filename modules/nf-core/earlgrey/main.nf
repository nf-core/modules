process EARLGREY {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/earlgrey:7.0.3--hd63eeec_0':
        'biocontainers/earlgrey:7.0.3--hd63eeec_0' }"

    input:
    tuple val(meta), path(genome)
    val(partition)

    output:
    tuple val(meta), path("${meta.id}"), emit: output_directory
    tuple val(meta), path("${meta.id}/*.fasta"), emit: fasta_files
    tuple val(meta), path("${meta.id}/*.gff"), emit: gff_files
    tuple val("${task.process}"), val('earlgrey'), eval('earlGrey -h | grep version | cut -d" " -f3'), topic: versions, emit: versions_earlgrey

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # unzip test
    gunzip *.gz

    # capture current working directory
    INPUT_DIR=\$PWD

    # first, change directory to the famdb library location
    cd /usr/local/share/RepeatMasker/Libraries/famdb/

    # download the partitions you require from Dfam 3.9
    curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.${partition}.h5.gz'

    # decompress Dfam 3.9 paritions
    gunzip *.gz

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

    earlGrey \\
        $args \\
        -g $genome \\
        -s ${prefix} \\
        -t $task.cpus \\
        -o ${prefix}
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
    echo $args

    mkdir ${prefix}
    """
}
