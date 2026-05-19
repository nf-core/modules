process RPBP_EXTRACTORFPROFILES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3a8aa95ce76934f6269b2d8cbdd3d57c13db029c704152975b2315e35b7a2154/data' :
        'community.wave.seqera.io/library/rpbp_star:247a8ae84a6babfb' }"

    input:
    tuple val(meta), path(bam), path(bai), path(periodic_offsets)
    path  orfs_genomic_bed
    path  exons_bed

    output:
    tuple val(meta), path("${prefix}.profiles.mtx.gz")             , emit: profiles
    tuple val(meta), path("${prefix}.periodic_lengths_offsets.tsv"), emit: lengths_offsets
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Periodic-length filter thresholds, space-separated:
    //   <min_metagene_profile_count> <min_metagene_bf_mean> <max_metagene_bf_var> <min_metagene_bf_likelihood>
    // Defaults mirror rpbp.defaults.metagene_options. Use "None" to disable a filter.
    def filter_args = (task.ext.args2 ?: '1000 5 None 0.5').tokenize(' ')
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_count      = filter_args[0]
    def min_bf_mean    = filter_args[1]
    def max_bf_var     = filter_args[2]
    def min_bf_lik     = filter_args[3]
    """
    # Replicates the length/offset filter from rpbp's
    # ribo_utils.utils.get_periodic_lengths_and_offsets so the upstream
    # `select-periodic-offsets` output can be consumed without going
    # through rpbp's filename-driven config plumbing. Same defaults as
    # rpbp.defaults.metagene_options.
    python - <<'PYTHON'
import pandas as pd, scipy.stats as st, math
df = pd.read_csv("${periodic_offsets}")
df = df[df["highest_peak_profile_sum"] >= ${min_count}]
df = df[df["highest_peak_bf_mean"]      >= ${min_bf_mean}]
max_bf_var = ${max_bf_var}
if max_bf_var is not None:
    df = df[df["highest_peak_bf_var"]   <= max_bf_var]
min_lik = ${min_bf_lik}
if min_lik is not None and len(df):
    s   = df["highest_peak_bf_var"].apply(math.sqrt)
    z   = (df["highest_peak_bf_mean"] - ${min_bf_mean}) / s
    lik = 1 - st.norm.cdf(-z)
    df  = df[lik > min_lik]
lengths = df["length"].astype(int).tolist()
offsets = df["highest_peak_offset"].astype(int).tolist()
with open("${prefix}.periodic_lengths_offsets.tsv", "w") as fh:
    fh.write("length\\toffset\\n")
    for l, o in zip(lengths, offsets):
        fh.write(f"{l}\\t{o}\\n")
print(" ".join(str(l) for l in lengths), file=open("lengths.txt","w"))
print(" ".join(str(o) for o in offsets), file=open("offsets.txt","w"))
PYTHON

    LENGTHS=\$(cat lengths.txt)
    OFFSETS=\$(cat offsets.txt)

    extract-orf-profiles \\
        ${bam} \\
        ${orfs_genomic_bed} \\
        ${exons_bed} \\
        ${prefix}.profiles.mtx.gz \\
        --lengths \$LENGTHS \\
        --offsets \$OFFSETS \\
        --num-cpus ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.profiles.mtx.gz
    touch ${prefix}.periodic_lengths_offsets.tsv
    """
}
