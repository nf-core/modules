- name: minimap2 align test_minimap2_align_single_end
  command: nextflow run tests/modules/minimap2/align -entry test_minimap2_align_single_end -c tests/config/nextflow.config
  tags:
    - minimap2
    - minimap2/align
  files:
    - path: output/minimap2/test.paf
      md5sum: 70e8cf299ee3ecd33e629d10c1f588ce
    - path: output/minimap2/versions.yml
      md5sum: cdb08919d584b7527241d191784b9da4

- name: y
  command: nextflow run tests/modules/minimap2/align -entry test_minimap2_align_paired_end -c tests/config/nextflow.config
  tags:
    - minimap2
    - minimap2/align
  files:
    - path: output/minimap2/test.paf
      md5sum: 5e7b55a26bf0ea3a2843423d3e0b9a28
    - path: output/minimap2/versions.yml
      md5sum: b748e1a56d7c638d3f1af495bb001c5d
