process GCTA_MAKEGRMPART {
    tag "part ${meta.part_gcta_job} of ${meta.nparts_gcta} (${meta.id})"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' }"

    input:
    tuple val(meta), path(mfile), path(bed_pgen), path(bim_pvar), path(fam_psam)
    tuple val(meta2), path(snp_group_file)

    output:
    tuple val(meta), path("${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.id"), path("${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.bin"), path("${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.N.bin"), emit: grm_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | grep 'version v' | tr -s ' ' | cut -d' ' -f3 | sed 's/^v//'"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def part_gcta_job = meta.part_gcta_job
    def nparts_gcta = meta.nparts_gcta
    def extract_cmd = snp_group_file ? "--extract ${snp_group_file}" : ''
    def extra_args = task.ext.args ?: ''
    def genotype_files = bed_pgen instanceof List ? bed_pgen : [bed_pgen]
    def genotype_extension = genotype_files[0].name.tokenize('.').last()
    def multi_file_flag = genotype_extension == 'pgen' ? '--mpfile' : '--mbfile'

    """
    set -euo pipefail

    gcta \\
        ${multi_file_flag} ${mfile} \\
        --make-grm-part ${nparts_gcta} ${part_gcta_job} \\
        ${extract_cmd} \\
        --maf 0.01 \\
        --thread-num ${task.cpus} \\
        --out ${meta.id} ${extra_args}
    """

    stub:
    """
    touch ${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.id
    touch ${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.bin
    touch ${meta.id}.part_${meta.nparts_gcta}_${meta.part_gcta_job}.grm.N.bin
    """
}
