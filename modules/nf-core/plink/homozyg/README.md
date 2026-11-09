# plink/homozyg

Runs of homozygosity (ROH) with PLINK 1.9 `--homozyg`.

This module only runs PLINK and writes its three output tables (`.hom`, `.hom.indiv`, `.hom.summary`) and log. The per-animal FROH, length-bin and plot steps are not part of this module; they are downstream analysis steps that are run in the pipeline that uses this output.

## Defaults

PLINK's own `--homozyg` defaults are used, so results match a plain PLINK run. The defaults are not tuned for a 60k goat array; pass explicit `--homozyg-*` options in `ext.args` to change them, and the values are recorded in the log.

## Downstream FROH

FROH per animal is total ROH length divided by the genome length (2,500 Mb by default in the pipeline). Animals must be keyed by FID + IID, because the same IID can occur under more than one population. Population values are averages of their animals' FROH, not sums.

## Validation

Tested on a goat SaudiGoatArdi input (144 animals, 59,727 SNPs, `--chr-set 31 --allow-extra-chr`): the `.hom`, `.hom.indiv` and `.hom.summary` outputs are identical to the recorded run (188 ROH segments).

## Container

The container is the biocontainers PLINK 1.90b6.21 image, the same as `plink/genome`.
