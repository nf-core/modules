process SEQKIT_SAMPLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit:2.13.0--205358a3675c7775'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'}"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("${prefix}*.${extension}"), emit: fastx
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^.*v//'"), emit: versions_seqkit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fastx_list = fastx instanceof List ? fastx : [fastx]
    def first_file = fastx_list[0].toString()
    extension = "fastq"
    if (first_file ==~ /.+\.(fasta|fa|fas|fna|fsa)(\.gz)?/) {
        extension = "fasta"
    }
    extension = first_file.endsWith('.gz') ? "${extension}.gz" : extension

    def cmds = []
    if (fastx_list.size() == 1) {
        def out_name = "${prefix}.${extension}"
        if (out_name == "${fastx_list[0]}") {
            error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
        }
        cmds << "seqkit sample --threads ${task.cpus} ${args} ${fastx_list[0]} -o ${out_name}"
    } else {
        fastx_list.eachWithIndex { f, i ->
            def out_name = "${prefix}_${i + 1}.${extension}"
            if (out_name == "${f}") {
                error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
            }
            cmds << "seqkit sample --threads ${task.cpus} ${args} ${f} -o ${out_name}"
        }
    }
    """
    ${cmds.join(' \\\n    && ')}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def fastx_list = fastx instanceof List ? fastx : [fastx]
    def first_file = fastx_list[0].toString()
    extension = "fastq"
    if (first_file ==~ /.+\.(fasta|fa|fas|fna|fsa)(\.gz)?/) {
        extension = "fasta"
    }
    extension = first_file.endsWith('.gz') ? "${extension}.gz" : extension
    def create_cmd = extension.endsWith('.gz') ? "echo '' | gzip >" : "touch"

    def stub_cmds = []
    if (fastx_list.size() == 1) {
        def out_name = "${prefix}.${extension}"
        if (out_name == "${fastx_list[0]}") {
            error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
        }
        stub_cmds << "${create_cmd} ${out_name}"
    } else {
        fastx_list.eachWithIndex { f, i ->
            def out_name = "${prefix}_${i + 1}.${extension}"
            if (out_name == "${f}") {
                error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
            }
            stub_cmds << "${create_cmd} ${out_name}"
        }
    }
    """
    ${stub_cmds.join('\n    ')}
    """
}
