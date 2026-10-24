// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_BED12CODONPOSITIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/7259cf5c00cd90b521cdda013955b69043542ad2b93f786c6ed53ecfc04a8422/data'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml:26e261815660ec9d'}"

    input:
    tuple val(meta), path(bed12)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml"                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''

    """
    echo ${args}
    """

    template 'bed12codonpositions.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

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
