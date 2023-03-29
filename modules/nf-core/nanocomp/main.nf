process NANOCOMP {
    label 'process_medium'

    conda "bioconda:nanocomp=1.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp:1.21.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanocomp:1.21.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(filelist)

    output:
    tuple val(meta), path("*.html"), emit: html_nanocomp_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    //determine input file type
    filetypes = []
    for (file in filelist){
        tokenized_filename = file.getName().tokenize('.')
        if (tokenized_filename.size() < 2){
            throw new java.lang.IndexOutOfBoundsException("Every input file to nanocomp has to have a file ending.")
        }
        
        first_namepart = true
        extension_found = false
        
        for (namepart in tokenized_filename){
            if (namepart == ""){
                continue
            }
            
            // prevent the file name to be seen as extension
            if (first_namepart == true){
                first_namepart = false
                continue
            }

            if (["fq","fastq"].contains(namepart)){
                filetypes.add("fastq")
                extension_found = true
                break
            } else if (["fasta", "fna", "ffn", "faa", "frn", "fa"].contains(namepart)) {
                filetypes.add("fasta")
                extension_found = true
                break
            } else if (namepart == "bam") {
                filetypes.add("fasta")
                extension_found = true
                break
            } else if (namepart == "txt") {
                filetypes.add("summary")
                extension_found = true
                break
            }
        }
        
        if (extension_found == false){
            throw new java.lang.IllegalArgumentException("There was no suitable filetype found for " + file.getName() + 
            ". NanoComp only accepts fasta (fasta, fna, ffn, faa, frn, fa), fastq (fastq, fq), bam and Nanopore sequencing summary (txt).")
        }
    }
    
    filetypes.unique()
    if (filetypes.size() < 1){
        throw new java.lang.IllegalArgumentException("There was no suitable filetype found in NanoComp input. Please use fasta, fastq, bam or Nanopore sequencing summary.")
    }
    if (filetypes.size() > 1){
        throw new java.lang.IllegalArgumentException("You gave different filetypes to NanoComp. Please use only *one* of fasta, fastq, bam or Nanopore sequencing summary.")
    }
    filetype = filetypes[0]

    """
    NanoComp \\
        --$filetype $filelist \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo \$(NanoComp --version 2>&1) | sed 's/^.*NanoComp //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
