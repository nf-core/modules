// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process METABINNER_METABINNER {
    tag "${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0':
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(kmer), path(depth)

    output:
    tuple val(meta), path("${prefix}.metabinner.tsv.gz"), emit: binning_results
    tuple val(meta), path("${prefix}.metabinner.log.gz"), emit: log
    tuple val("${task.process}"), val('metabinner'), val('1.4.4'), emit: versions_metabinner, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def (assembly_file, kmer_file, depth_file) = [assembly, kmer, depth].collect { file -> file.name.replaceFirst(/\.gz$/, '') }
    def coverage_profile = "${prefix}.metabinner_coverage.tsv"
    """
    for file in ${assembly} ${kmer} ${depth}; do
        if [[ "\$file" == *.gz ]]; then
            gzip -cd "\$file" > "\${file%.gz}"
        fi
    done

    # MetaBinner bins every contig in the coverage profile, so restrict it to the k-mer profile's contigs
    # (and thereby its length threshold); keep only per-sample mean depths, whatever the jgi options
    touch ${prefix}.warnings.txt
    awk -F'\\t' -v OFS='\\t' '
            NR == FNR { if (FNR > 1) in_kmer[\$1] = 1; next }
            FNR == 1 { for (i = 1; i <= NF; i++) if (i == 1 || (\$i != "contigLen" && \$i != "totalAvgDepth" && \$i !~ /-var\$/)) keep[++n] = i }
            FNR == 1 || \$1 in in_kmer { line = \$keep[1]; for (j = 2; j <= n; j++) line = line OFS \$keep[j]; print line; next }
            { dropped++ }
            END { if (dropped) printf "WARNING: %d contig(s) in ${depth} are absent from ${kmer} and were not binned (e.g. below the k-mer profile length threshold)\\n", dropped > "${prefix}.warnings.txt" }
        ' FS=',' ${kmer_file} FS='\\t' ${depth_file} > ${coverage_profile}
    cat ${prefix}.warnings.txt >&2

    if [ \$(wc -l < ${coverage_profile}) -le 1 ]; then
        echo "ERROR: no contigs shared between ${depth} and ${kmer}. Do both come from the same assembly?" >&2
        exit 1
    fi

    metabinner_path=\$(dirname \$(which run_metabinner.sh))

    # run_metabinner.sh requires absolute paths
    run_metabinner.sh \\
        -a \$PWD/${assembly_file} \\
        -k \$PWD/${kmer_file} \\
        -d \$PWD/${coverage_profile} \\
        -o \$PWD/${prefix} \\
        -t ${task.cpus} \\
        -p \$metabinner_path \\
        ${args}

    # MetaBinner output order is non-deterministic
    LC_ALL=C sort -k1,1 ${prefix}/metabinner_res/metabinner_result.tsv | gzip -cn > ${prefix}.metabinner.tsv.gz
    cat ${prefix}.warnings.txt ${prefix}/metabinner_res/result.log > ${prefix}.metabinner.log
    gzip -n ${prefix}.metabinner.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    echo $args
    echo "" | gzip > ${prefix}.metabinner.tsv.gz
    echo "" | gzip > ${prefix}.metabinner.log.gz
    """
}
