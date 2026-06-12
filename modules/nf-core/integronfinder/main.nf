process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_low'

    containerOptions workflow.containerEngine in ['singularity', 'apptainer'] ?
        '--env MPLCONFIGDIR=/tmp/mplconfigdir' :
        '-e MPLCONFIGDIR=/tmp/mplconfigdir'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/integron_finder:2.0.5--pyhdfd78af_0':
        'quay.io/biocontainers/integron_finder:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*/*.gbk")              , emit: gbk, optional: true
    tuple val(meta), path("*/*.integrons")        , emit: integrons
    tuple val(meta), path("*/*.summary")          , emit: summary
    tuple val(meta), path("*/integron_finder.out"), emit: out
    tuple val("${task.process}"), val('integronfinder'), eval("integron_finder --version 2>&1 | sed '1!d;s/^integron_finder version //;s/ .*//'"), emit: versions_integronfinder, topic: versions


    script:
    def args                = task.ext.args ?: ''
    def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name          = fasta.getName().replace(".gz", "")

    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    integron_finder \\
        $args \\
        --cpu $task.cpus \\
        $fasta_name
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gbk_create = args.contains("--gbk") ? "Results_Integron_Finder_${prefix}/${prefix}.gbk" : ""

    """
    mkdir -p Results_Integron_Finder_${prefix}
    ${gbk_create}
    touch "Results_Integron_Finder_${prefix}/${prefix}.integrons"
    touch "Results_Integron_Finder_${prefix}/${prefix}.summary"
    touch "Results_Integron_Finder_${prefix}/integron_finder.out"
    """
}
