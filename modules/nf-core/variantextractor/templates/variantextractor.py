#!/usr/bin/env python3

import pysam
import variant_extractor
from variant_extractor import VariantExtractor

vcf_in = pysam.VariantFile("${vcf}")
header = str(vcf_in.header).rstrip("\\n")
vcf_in.close()

pass_only_str = "${pass_only}"
ensure_pairs_str = "${ensure_pairs}"

extractor = VariantExtractor(
    "${vcf}",
    pass_only=pass_only_str == "True",
    ensure_pairs=ensure_pairs_str == "True",
)

with open("${prefix}.extracted.vcf", "w") as out:
    out.write(header + "\\n")
    for variant_record in extractor:
        out.write(str(variant_record) + "\\n")

extractor.close()

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f"    variant-extractor: {variant_extractor.__version__}\\n")
