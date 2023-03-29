- name: cellpose test_cellpose
  command: nextflow run ./tests/modules/nf-core/cellpose -entry test_cellpose -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/cellpose/nextflow.config
  tags:
    - cellpose
  files:
    - path: output/cellpose/cycif_tonsil_small.ome_cp_masks.tif
      md5sum: 59af5f94e1bf55a88b50a840c6cfb244
    - path: output/cellpose/versions.yml
