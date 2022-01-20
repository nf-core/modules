process MAXQUANT_LFQ {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::maxquant=2.0.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/maxquant:2.0.1.0--py39hdfd78af_2"
    } else {
	container "wombatp/maxquant-pipeline:v0.13"
    }

    input:
    tuple val(meta), path(fasta), path(paramfile)
    path raw

    output:
    tuple val(meta), path("*.txt"), emit: maxquant_txt
    path "versions.yml"          , emit: version

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"  
    
    """
    export PATH=/usr/local/lib/dotnet:/usr/local/lib/dotnet/tools:/opt/conda/envs/nf-core-maxquant/bin:/opt/conda/envs/nf-core-maxquant/lib/dotnet/tools:/opt/conda/envs/nf-core-maxquant/lib/dotnet:$PATH	
    echo \"maxquant: \"\$(maxquant --version 2>&1 > /dev/null | cut -f2 -d\" \") > versions.yml
    sed \"s_<numThreads>.*_<numThreads>$task.cpus</numThreads>_\" ${paramfile} > mqpar_changed.xml
    sed -i \"s|PLACEHOLDER|\$PWD/|g\" mqpar_changed.xml
    mkdir temp
    maxquant mqpar_changed.xml
    mv combined/txt/*.txt .
    """
}
