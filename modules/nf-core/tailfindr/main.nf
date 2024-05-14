process TAILFINDR {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ont-fast5-api=0.4.1 bioconda::r-tailfindr=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f24f1399a77784f913670cbb36a0f17b78e0631b:80e40d512cd5a71665e3e00e8d0ad1462fc58f76-0':
        'biocontainers/mulled-v2-f24f1399a77784f913670cbb36a0f17b78e0631b:80e40d512cd5a71665e3e00e8d0ad1462fc58f76-0' }"

    input:
    tuple val(meta), path(fast5)

    output:
    tuple val(meta), path("*.csv.gz"), emit: csv_gz
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? ", ${task.ext.args}": ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R --vanilla --slave -e "library(tailfindr);
    find_tails(fast5_dir = './' ,
    save_dir = './' ${args},
    csv_filename = \'${meta.id}.csv\',
    num_cores = ${task.cpus})";
    gzip ${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tailfindr: \$(Rscript -e "cat(paste(packageVersion('tailfindr'), collapse='.'))")
        ont-fast5-api: \$(pip show ont-fast5-api | grep Version | awk '{print \$2}')
    END_VERSIONS
    """
}
