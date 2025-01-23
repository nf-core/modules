import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { BWA_INDEX  } from '/home/ec2-user/git/modules/modules/nf-core/parabricks/fq2bam/tests/../../../bwa/index/main.nf'

include { BWA_INDEX  as BWA_INDEX_PE } from '/home/ec2-user/git/modules/modules/nf-core/parabricks/fq2bam/tests/../../../bwa/index/main.nf'


// include test process
include { PARABRICKS_FQ2BAM } from '/home/ec2-user/git/modules/modules/nf-core/parabricks/fq2bam/tests/../main.nf'

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
                    [ id:'test' ], // meta map
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true)
                ])
                
        BWA_INDEX(*input)
    }
    
    {
        def input = []
        
                input[0] = Channel.of([
                    [ id:'test' ], // meta map
                    file(params.modules_testdata_base_path + 'genomics/sarscov2/genome/genome.fasta', checkIfExists: true)
                ])
                
        BWA_INDEX_PE(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = Channel.of([
				    [ 'id':'test', 'read_group':'rg', 'sample':'sm', 'platform':'pl', 'single_end':true ],
                    [ [ 'read_group':'rg', 'sample':'sm', 'platform':'pl' ] ],
				    [
				        file('https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz', checkIfExists: true)
				    ],
                    []
				])
                input[1] = Channel.of([
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa', checkIfExists: true),
                    file('https://github.com/nf-core/test-datasets/raw/methylseq/reference/genome.fa.fai', checkIfExists: true)
                ])
				input[2] = BWA_INDEX.out.index
                input[3] = []
                input[4] = []
                
    //----

    //run process
    PARABRICKS_FQ2BAM(*input)

    if (PARABRICKS_FQ2BAM.output){

        // consumes all named output channels and stores items in a json file
        for (def name in PARABRICKS_FQ2BAM.out.getNames()) {
            serializeChannel(name, PARABRICKS_FQ2BAM.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = PARABRICKS_FQ2BAM.out as Object[]
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
