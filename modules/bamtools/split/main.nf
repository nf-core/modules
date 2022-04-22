process BAMTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bamtools=2.5.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hd03093a_0' :
        'quay.io/biocontainers/bamtools:2.5.2--hd03093a_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_list = bam.collect{"-in $it"}.join(' ')
    def stub_cmd = !args.contains("-stub") : " -stub ${prefix}" : ""

    if (size(bam) > 1){
        def bamtools_merge_cmd = "bamtools merge ${input_list} |"
    } else {
        def split_input = input_list
    }


    """
    ${bamtools_merge_cmd} \\
    bamtools split \\
        $split_input \\
        $stub_cmd \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
