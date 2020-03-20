nextflow.preview.dsl=2
params.genome = ''

process HISAT2 {
    // depending on the genome used one might want/need to adjust the memory settings.
    // For the E. coli test data this is probably not required
    // label 'bigMem'
    // label 'multiCore'

    input:
        tuple val(name), path(reads)
        val (outdir)
        val (hisat2_args)
        val (verbose)

    output:
        path "*bam",       emit: bam
        path "*stats.txt", emit: stats 

    publishDir "$outdir/hisat2",
        mode: "copy", overwrite: true

    script:
    
        if (verbose){
            println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
        }
    
        cores = 4
        readString = ""
        hisat_options = hisat2_args

        // Options we add are
        hisat_options = hisat_options + " --no-unal --no-softclip "

        if (reads instanceof List) {
            readString = "-1 "+reads[0]+" -2 "+reads[1]
            hisat_options = hisat_options + " --no-mixed --no-discordant"
        }
        else {
            readString = "-U "+reads
        }
        index = params.genome["hisat2"]
        
        splices = ''
        if (params.genome.containsKey("hisat2_splices")){
            splices = " --known-splicesite-infile " + params.genome["hisat2_splices"]
        }
        else{
            println ("No key 'hisat2_splices' was supplied. Skipping...")
        }
        hisat_name = name + "_" + params.genome["name"]

        """
        module load hisat2
        module load samtools
        hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
        """

}