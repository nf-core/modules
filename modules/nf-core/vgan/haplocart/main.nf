def VERSION = '1.0.1'

process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vgan=1.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:1.0.1--h9ee0642_0':
        'quay.io/biocontainers/vgan' }"

    input:
    tuple val(meta), path(reads), path(reads2)

    output:
    tuple val(meta), path("*.txt"), path("*.posterior.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    mkdir hc_files
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.gg"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.og"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.xg"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.giraffe.gbz"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.dist"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph_paths"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/path_supports"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/graph.gbwt"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/children.txt"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/parents.txt"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/parsed_pangenome_mapping"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/mappability.tsv"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/k17_w18.min"
    wget -nc --no-parent -P hc_files "ftp://ftp.healthtech.dtu.dk/public/haplocart/hcfiles/k31_w11.min"


    if [[ "${meta.format}" == "fastq" ]] && ${meta.single_end};
    then 
        vgan haplocart $args -t $task.cpus -fq1 $reads -o ${prefix}.txt --hc-files hc_files -pf ${prefix}.posterior.txt;
    elif [[ "${meta.format}" == "fastq" ]] && [[ "$reads2" != "" ]];
    then
        vgan haplocart $args -t $task.cpus -fq1 $reads -fq2 $reads2 -o ${prefix}.txt --hc-files hc_files -pf ${prefix}.posterior.txt;
    elif [[ "${meta.format}" == "fastq" ]];
    then
        vgan haplocart $args -t $task.cpus -fq1 $reads -i -o ${prefix}.txt --hc-files hc_files -pf ${prefix}.posterior.txt;
    fi

    echo $VERSION >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.posterior.txt

    echo $VERSION >> versions.yml
    """

}
