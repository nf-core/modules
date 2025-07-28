import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { MINIMAP2_ALIGN  } from '/workspace/modules/modules/nf-core/comebin/runcomebin/tests/../../../minimap2/align/main.nf'


// include test process
include { COMEBIN_RUNCOMEBIN } from '/workspace/modules/modules/nf-core/comebin/runcomebin/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    
    {
        def input = []
        
                    input[0] = Channel.of([
                        [id: "test"],
                        [
                            file("https://tolit.cog.sanger.ac.uk/test-data/Zymo_D6311_Metagenome/genomic_data/pacbio/fasta/subsampled.pacbio.part_001.fa.gz", checkIfExists: true),
                            file("https://tolit.cog.sanger.ac.uk/test-data/Zymo_D6311_Metagenome/genomic_data/pacbio/fasta/subsampled.pacbio.part_002.fa.gz", checkIfExists: true)
                        ]
                    ])
                    input[1] = [
                        [id: "ref"],
                        file(
                            "https://tolit.cog.sanger.ac.uk/test-data/Zymo_D6311_Metagenome/assembly/xyTesTing1_metamdbg.contigs.fasta.gz",
                            checkIfExists: true
                        )
                    ]
                    input[2] = true
                    input[3] = "csi"
                    input[4] = false
                    input[5] = false
                    
        MINIMAP2_ALIGN(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = Channel.of([
                        [id: "ref"],
                        file(
                            "https://tolit.cog.sanger.ac.uk/test-data/Zymo_D6311_Metagenome/assembly/xyTesTing1_metamdbg.contigs.fasta.gz",
                            checkIfExists: true
                        )
                    ]).combine(MINIMAP2_ALIGN.out.bam)
                    .map { ref_meta, asm, reads_meta, bam -> [reads_meta, asm, bam] }
                
    //----

    //run process
    COMEBIN_RUNCOMEBIN(*input)

    if (COMEBIN_RUNCOMEBIN.output){

        // consumes all named output channels and stores items in a json file
        for (def name in COMEBIN_RUNCOMEBIN.out.getNames()) {
            serializeChannel(name, COMEBIN_RUNCOMEBIN.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = COMEBIN_RUNCOMEBIN.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
