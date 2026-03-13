cd /workspaces/modules

#nf-test test modules/nf-core/cellranger/count  --profile docker --update-snapshot
#nf-test test modules/nf-core/cellranger/mkgtf  --profile docker --update-snapshot
#nf-test test modules/nf-core/cellranger/mkref  --profile docker --update-snapshot
#nf-test test modules/nf-core/cellranger/mkvdjref  --profile docker --update-snapshot
#nf-test test modules/nf-core/cellranger/multi  --profile docker --update-snapshot
nf-test test modules/nf-core/cellranger/vdj  --profile docker --update-snapshot
