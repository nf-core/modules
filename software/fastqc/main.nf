nextflow.preview.dsl = 2
def MODULE = "fastqc"
params.publish_dir = MODULE
params.publish_results = "default"

process FASTQC {
    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path ("*fastqc*"), emit: all
        path "*.zip", emit: report // e.g. for MultiQC later
        path "*.version.txt", emit: version

    container "docker.pkg.github.com/nf-core/$MODULE"
    conda "${moduleDir}/environment.yml"

    publishDir "${params.out_dir}/${params.publish_dir}/$name",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                        if(params.publish_results == "none") null
                        else filename }

    script:
        """
        fastqc ${params.fastqc_args} -t ${task.cpus} $reads
        fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
        """
}
