process SPIKEINTERFACE_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spikeinterface:0.101.2--pyh0610db2_0' :
        'biocontainers/spikeinterface:0.101.2--pyh0610db2_0' }"

    input:
    tuple val(meta), path(recording)
    val   sorter

    output:
    tuple val(meta), path("${prefix}/sorting"),                    emit: sorting
    tuple val(meta), path("${prefix}/sorting/spike_times.npy"),    emit: spike_times
    tuple val(meta), path("${prefix}/sorting/spike_clusters.npy"), emit: spike_clusters
    path  "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    python - <<'PY'
import os, numpy as np
import spikeinterface.extractors as se
import spikeinterface.sorters as ss

recording = se.read_binary_folder("${recording}") if os.path.isdir("${recording}") \\
    else se.read_binary("${recording}", sampling_frequency=30000.0,
                        dtype="int16", num_channels=1)

sorting = ss.run_sorter(
    sorter_name="${sorter}",
    recording=recording,
    folder="${prefix}/sorting",
    remove_existing_folder=True,
)

spikes = sorting.get_all_spike_trains()[0]
np.save("${prefix}/sorting/spike_times.npy",    spikes[0])
np.save("${prefix}/sorting/spike_clusters.npy", spikes[1])
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spikeinterface: 0.101.2
        spike_sorter: ${sorter}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/sorting
    touch ${prefix}/sorting/spike_times.npy
    touch ${prefix}/sorting/spike_clusters.npy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spikeinterface: 0.101.2
        spike_sorter: ${sorter}
    END_VERSIONS
    """
}
