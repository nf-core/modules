process GCTA_MAKEGRM_HELPER {
    tag "${meta.id}"
    label "process_medium"
    conda "${projectDir}/modules/nf-core/gcta/bivariateremlldms/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9' :
        'community.wave.seqera.io/library/gcta:1.94.1--9bc35dc424fcf6e9'}"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(extract_file)

    output:
    tuple val(meta), path("${meta.id}.grm.id"), path("${meta.id}.grm.bin"), path("${meta.id}.grm.N.bin"), emit: grm_files

    script:
    def extract_param = extract_file ? "--extract ${extract_file}" : ""
    def bfile_prefix = bed.baseName
    """
    set -euo pipefail

    gcta \\
        --bfile "${bfile_prefix}" \\
        ${extract_param} \\
        --make-grm \\
        --out "${meta.id}" \\
        --thread-num ${task.cpus}
    """

    stub:
    """
    touch "${meta.id}.grm.id"
    touch "${meta.id}.grm.bin"
    touch "${meta.id}.grm.N.bin"
    """
}
