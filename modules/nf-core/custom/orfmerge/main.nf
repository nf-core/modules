process CUSTOM_ORFMERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1019bd22c111267bcb670fdb128829776f0ca6adfa7b0e2d126f91577d08e3/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:75514f9f977be607' }"

    input:
    tuple val(meta), path(bed12s, arity: '1..*', stageAs: 'beds/*'), path(tsvs, arity: '1..*', stageAs: 'tsvs/*')

    output:
    tuple val(meta), path("${prefix}.catalogue.bed12")    , emit: bed12
    tuple val(meta), path("${prefix}.catalogue.tsv")      , emit: catalogue_tsv
    tuple val(meta), path("${prefix}.orf_to_gene.tsv")    , emit: orf_to_gene_tsv
    tuple val(meta), path("${prefix}.catalogue.mqc.tsv")  , emit: multiqc
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''
    template 'orfmerge.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.catalogue.bed12
    touch ${prefix}.catalogue.tsv
    touch ${prefix}.orf_to_gene.tsv
    touch ${prefix}.catalogue.mqc.tsv

    python - <<END
import platform
import pandas
import yaml

with open("versions.yml", "w") as fh:
    yaml.safe_dump(
        {
            "${task.process}": {
                "python": platform.python_version(),
                "pandas": pandas.__version__,
            }
        },
        fh,
        default_flow_style=False,
        sort_keys=False,
    )
END
    """
}
