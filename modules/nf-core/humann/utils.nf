def getProcessName(task_process) {
    return task_process.tokenize(':')[-1]
}

def getContainer(name)  {
    return [
	'HUMANN3': 'ghcr.io/vdblab/biobakery-profiler:4.0.5--3.6.1_smaller-pt2',
	'HUMANN4': 'ghcr.io/vdblab/biobakery-profiler:4.0.6--4.0.0.alpha.1-final_smaller-pt2'
    ][name]
}
def getConda(name) {
    return [
	'HUMANN3': 'bioconda::humann=3.6.1',
	'HUMANN4': 'bioconda::humann=4.0.0.alpha.1-final'
    ][name]
}
def getExt(name) {
    return [
	'HUMANN3': '*.ffn.gz',
	'HUMANN4': '*.fna.gz'
    ][name]
}



// These mock processes are here to ensure that both containers get downloaded by nextflow inspect/download.  Without this, the command doesn't evaluate the aliases which we are using to select the proper containers
process testhumann3{
    container  'ghcr.io/vdblab/biobakery-profiler:4.0.5--3.6.1_smaller-pt2'

    output:
    val "x3"
    script: """echo x3 """
}
process testhumann4{
    container  'ghcr.io/vdblab/biobakery-profiler:4.0.6--4.0.0.alpha.1-final_smaller-pt2'
    output:
    val "x4"
    script: """echo x4 """

}
