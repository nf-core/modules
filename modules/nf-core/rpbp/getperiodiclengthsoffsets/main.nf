process RPBP_GETPERIODICLENGTHSOFFSETS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/146c3f15abf184a5ec13531d2a040ba7b9235c1091723aa37c7a119817411367/data' :
        'community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b' }"

    input:
    tuple val(meta), path(periodic_offsets)

    output:
    tuple val(meta), path("${prefix}.periodic_lengths_offsets.tsv"), emit: lengths_offsets
    tuple val("${task.process}"), val('rpbp'), eval('python -c "import rpbp; print(rpbp.__version__)"'), emit: versions_rpbp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Periodic-length filter thresholds, space-separated:
    //   <min_metagene_profile_count> <min_metagene_bf_mean> <max_metagene_bf_var> <min_metagene_bf_likelihood>
    // Defaults mirror rpbp.defaults.metagene_options. Any token may be "None"
    // to defer to upstream's per-slot default-filter handling.
    def filter_args = (task.ext.args ?: '1000 5 None 0.5').tokenize(' ')
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_count   = filter_args[0]
    def min_bf_mean = filter_args[1]
    def max_bf_var  = filter_args[2]
    def min_bf_lik  = filter_args[3]
    """
    mkdir -p rpbp_work/metagene-profiles
    cp ${periodic_offsets} rpbp_work/metagene-profiles/${prefix}-unique.periodic-offsets.csv.gz

    python <<'PYTHON'
import pandas as pd
from rpbp.ribo_utils.utils import get_periodic_lengths_and_offsets

config = dict(
    riboseq_data="rpbp_work",
    min_metagene_profile_count=${min_count},
    min_metagene_bf_mean=${min_bf_mean},
    max_metagene_bf_var=${max_bf_var},
    min_metagene_bf_likelihood=${min_bf_lik},
)
lengths, offsets = get_periodic_lengths_and_offsets(config, "${prefix}", is_unique=True)
pd.DataFrame({"length": lengths, "offset": offsets}).to_csv(
    "${prefix}.periodic_lengths_offsets.tsv", sep="\\t", index=False)
PYTHON
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.periodic_lengths_offsets.tsv
    """
}
