import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { DEEPTOOLS_BIGWIGCOMPARE } from '/home/nadiunix/modules/modules/nf-core/deeptools/bigwigcompare/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                def bigwig1 = file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bigwig/test_S2.RPKM.bw', checkIfExists: true)
                def bigwig2 = file(params.modules_testdata_base_path + 'genomics/homo_sapiens/illumina/bigwig/test_S3.RPKM.bw', checkIfExists: true).copyTo('test2.bigwig')

                input[0] = [
                    [ id:'test' ],
                    bigwig1, bigwig2
                ]
                input[1] = []
                
    //----

    //run process
    DEEPTOOLS_BIGWIGCOMPARE(*input)

    if (DEEPTOOLS_BIGWIGCOMPARE.output){

        // consumes all named output channels and stores items in a json file
        for (def name in DEEPTOOLS_BIGWIGCOMPARE.out.getNames()) {
            serializeChannel(name, DEEPTOOLS_BIGWIGCOMPARE.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = DEEPTOOLS_BIGWIGCOMPARE.out as Object[]
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
