process GCTA_CALCULATELDSCORES {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}_gcta_ld.score.ld"), emit: ld_scores
    tuple val(meta), path("${meta.id}_snp_group1.txt"), path("${meta.id}_snp_group2.txt"), path("${meta.id}_snp_group3.txt"), path("${meta.id}_snp_group4.txt"), emit: snp_group_files
    tuple val("${task.process}"), val("gcta"), eval("gcta --version 2>&1 | head -n 1"), emit: versions_gcta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def extra_args = task.ext.args ?: ''

    """
    set -euo pipefail

    gcta \\
        --bfile ${meta.id} \\
        --ld-score-region 200 \\
        --out ${meta.id}_gcta_ld \\
        --thread-num ${task.cpus} ${extra_args}

    ld_file="${meta.id}_gcta_ld.score.ld"
    sorted_file="ldscore.sorted.tsv"

    awk 'NR > 1 { print \$1 "\\t" \$8 }' "${meta.id}_gcta_ld.score.ld" | sort -k2,2n > "\${sorted_file}"

    count=\$(wc -l < "\${sorted_file}")
    q1_idx=\$(( (count + 3) / 4 ))
    q2_idx=\$(( (count + 1) / 2 ))
    q3_idx=\$(( (3 * count + 1) / 4 ))

    q1=\$(awk -v idx="\${q1_idx}" 'NR == idx { print \$2 }' "\${sorted_file}")
    q2=\$(awk -v idx="\${q2_idx}" 'NR == idx { print \$2 }' "\${sorted_file}")
    q3=\$(awk -v idx="\${q3_idx}" 'NR == idx { print \$2 }' "\${sorted_file}")

    awk -v q1="\${q1}" -v q2="\${q2}" -v q3="\${q3}" -v prefix="${meta.id}" '
    NR > 1 {
        if (\$8 <= q1) {
            print \$1 >> prefix "_snp_group1.txt"
        } else if (\$8 <= q2) {
            print \$1 >> prefix "_snp_group2.txt"
        } else if (\$8 <= q3) {
            print \$1 >> prefix "_snp_group3.txt"
        } else {
            print \$1 >> prefix "_snp_group4.txt"
        }
    }
    ' "\${ld_file}"
    """

    stub:
    """
    touch ${meta.id}_gcta_ld.score.ld
    touch ${meta.id}_snp_group1.txt
    touch ${meta.id}_snp_group2.txt
    touch ${meta.id}_snp_group3.txt
    touch ${meta.id}_snp_group4.txt
    """
}
