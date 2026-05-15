process GCTA_MAKEGRM {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46b0d05f0daa47561d87d2a9cac5e51edc2c78e26f1bbab439c688386241a274/data'
        : 'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(mfile), path(bed_pgen), path(bim_pvar), path(fam_psam)

    output:
    tuple val(meta), path("*.grm.*"), emit: grm_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version | sed -En 's/^[*] version v([0-9.]*).*/\\1/p'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genotype_files = bed_pgen instanceof List ? bed_pgen : [bed_pgen]
    def genotype_extension = genotype_files[0].name.tokenize('.').last()
    def multi_file_flag = genotype_extension == 'pgen' ? '--mpfile' : '--mbfile'

    """

    gcta \\
        ${multi_file_flag} ${mfile} \\
        --make-grm \\
        --thread-num ${task.cpus} \\
        --out ${prefix} ${extra_args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.grm.id
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.N.bin
    """
}
