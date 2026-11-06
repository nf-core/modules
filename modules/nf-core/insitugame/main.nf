// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process INSITUGAME {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/aa/aa25d488809b83440425f792955c529c86887805184f54ec639323723b8bf930/data'
        : 'community.wave.seqera.io/library/insitugame:1.0.7--8a968e208a8a3d71'}"

    input:
    tuple val(meta), path(spatial_data)

    output:
    tuple val(meta), path("${prefix}_results.json"), emit: results
    tuple val(meta), path("${prefix}_tables.csv"), emit: tables
    tuple val("${task.process}"), val('insitugame'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('insitugame'))\""), topic: versions, emit: versions_insitugame

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c 'import csv, gzip, pathlib, sys; p=pathlib.Path("${spatial_data}"); f=(gzip.open if p.suffix == ".gz" else open)(p, "rt", newline=""); first=f.readline(); f.seek(0); h=next(csv.reader(f, delimiter="\\t" if "\\t" in first else ",")); required={"T_S":{"T_S","Sensitive"},"T_R":{"T_R","Resistant"},"S_CAF":{"S_CAF","Stroma"},"FITNESS":{"FITNESS","Fitness"},"HYPOXIA":{"HYPOXIA","Hypoxia"},"VASCULAR":{"VASCULAR","Vascular"}}; missing=[key for key, names in required.items() if not names.intersection(h)]; sys.exit("Missing required marker columns: " + ", ".join(missing) if missing else 0)'
    insitugame \\
        ${args} \\
        --input "${spatial_data}" \\
        --output-prefix "${prefix}"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.json
    touch ${prefix}_tables.csv
    """
}
