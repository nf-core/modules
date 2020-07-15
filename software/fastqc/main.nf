def MODULE = "fastqc"
params.publish_dir = MODULE
params.publish_results = "default"

process FASTQC {
    publishDir "${params.out_dir}/${params.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (params.publish_results == "none") null
                    else filename }

    container "docker.pkg.github.com/nf-core/$MODULE"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(name), val(single_end), path(reads)

    output:
    tuple val(name), val(single_end), path("*.html"), emit: html
    tuple val(name), val(single_end), path("*.zip"), emit: zip
    path "*.version.txt", emit: version

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    if (single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        fastqc ${params.fastqc_args} --threads $task.cpus ${name}.fastq.gz
        fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        fastqc ${params.fastqc_args} --threads $task.cpus ${name}_1.fastq.gz ${name}_2.fastq.gz
        fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
        """
    }
}
