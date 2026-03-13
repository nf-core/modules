process VARIANTEXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/variant-extractor:5.1.0--pyh106432d_0':
        'biocontainers/variant-extractor:5.1.0--pyh106432d_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('variant-extractor'), eval("python3 -c \"import variant_extractor; print(variant_extractor.__version__)\""), topic: versions, emit: versions_variantextractor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def pass_only    = args.contains('--pass-only')       ? 'True'  : 'False'
    def ensure_pairs = args.contains('--no-ensure-pairs') ? 'False' : 'True'
    """
    python3 << 'PYEOF'
import pysam
from variant_extractor import VariantExtractor

vcf_in   = pysam.VariantFile("${vcf}")
header   = str(vcf_in.header).rstrip("\\n")
vcf_in.close()

extractor = VariantExtractor(
    "${vcf}",
    pass_only=${pass_only},
    ensure_pairs=${ensure_pairs}
)

with open("${prefix}.extracted.vcf", "w") as out:
    out.write(header + "\\n")
    for variant_record in extractor:
        out.write(str(variant_record) + "\\n")

extractor.close()
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.extracted.vcf
    """
}
