process HTSLIB_REBGZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33a1f2c7f36ec58339e41cbea096d121f606918778a91cfbef944b40ba7ce48b/data'
        : 'community.wave.seqera.io/library/htslib_xz:49c8c84af5c4b3b9'}"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("out/${outfile}"), emit: bgzipped
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1! d; s/bgzip (htslib) //'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def in_ext = "${infile.extension}"
    def zip_extensions = ["xz", "bz2", "gz"]
    def prefix = task.ext.prefix ?: zip_extensions.contains(in_ext.toString()) ? "${infile.baseName}" : "${infile.name}"
    in_cmd = in_ext == "gz" ? "zcat" : ( in_ext == "xz" ? "xzcat" : ( in_ext == "bz2" ? "bzcat" : "cat" ) )
    outfile = "${prefix}.gz"
    """
    mkdir -p out
    ${in_cmd} ${infile} | bgzip -c ${args} -@ ${task.cpus} > "out/${outfile}"
    """

    stub:
    def args = task.ext.args ?: ''
    in_ext = "${infile.extension}"
    in_cmd = in_ext == "gz" ? "zcat" : ( in_ext == "xz" ? "xzcat" : ( in_ext == "bz2" ? "bzcat" : "cat" ) )
    outname = ["xz", "bz2", "gz"].contains(in_ext) ? "${infile.baseName}" : "${infile.name}"
    outfile = task.ext.prefix ? "${task.ext.prefix}.gz" : "${outname}.gz"
    """
    echo "$in_cmd" | bgzip -c ${args} -@ ${task.cpus} > "out/${outfile}"
    """
}
