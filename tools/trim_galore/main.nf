nextflow.preview.dsl=2

params.singlecell = ''
params.rrbs       = ''
params.pbat       = ''
params.single_end  = false

params.trim_nextseq = 0

params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0


process TRIM_GALORE {	
    
    // container 'quay.io/biocontainers/trim-galore:0.6.5--0' // maybe later
    // tag "$sample_id"

	input:
	    tuple val (name), path (reads)
		val (outdir)
		val (trim_galore_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fq.gz"),             emit: reads
		path "*trimming_report.txt", optional: true,  emit: report
		
        // Trimming reports are not generated for e.g. --hardtrim5, --clock etc
        // saveAs: {filename ->
        //   else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
        //   else filename
        // }

	publishDir "${outdir}/trim_galore",
		mode: "copy", overwrite: true

    script:
		if (verbose){
			println ("[MODULE] TRIM GALORE ARGS: " + trim_galore_args)
		}
		
        trim_galore_args += " --gzip "           // we like small files

		pairedString = 0
		if (reads instanceof List) {
			pairedString = 1
            trim_galore_args += " --paired "
		}
		
        if (params.clip_r1 > 0){
            trim_galore_args += " --clip_r1 ${params.clip_r1} "
        }
        if (params.clip_r2 > 0){
            trim_galore_args += " --clip_r2 ${params.clip_r2} "
        }
        if (params.three_prime_clip_r1> 0){
            trim_galore_args += " --three_prime_clip_r1 ${params.three_prime_clip_r1} "
        }
        if (params.three_prime_clip_r2 > 0){
            trim_galore_args += " --three_prime_clip_r2 ${params.three_prime_clip_r2} "
        }
        
        if (params.trim_nextseq > 0){
            trim_galore_args += " --nextseq ${params.trim_nextseq} "
        }   
    
        
        // Pre-set parameters for certain bisulfite-seq applications
        if (params.singlecell){
            trim_galore_args += " --clip_r1 6 "
            if (pairedString == 1){
                trim_galore_args += " --clip_r2 6 "
            }
        }
        if (params.rrbs){
            trim_galore_args += " --rrbs "
        } 
        if  (params.pbat){
            trim_galore_args += " --clip_r1 $params.pbat "
            if (pairedString == 1){
                trim_galore_args += " --clip_r2 $params.pbat "
            }
        }

		"""
		module load trim_galore
		trim_galore $trim_galore_args $reads
		"""

}



 


  

