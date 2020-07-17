#!/usr/bin/env nextflow
import java.security.MessageDigest
nextflow.preview.dsl = 2

params.out_dir = "test_output"
params.fastqc_args = ''
params.publish_dir_mode = "copy"

include { FASTQC } from '../main.nf'


private static String getFileChecksum(MessageDigest digest, File file) throws IOException
{
    // https://howtodoinjava.com/java/io/how-to-generate-sha-or-md5-file-checksum-hash-in-java/
    //Get file input stream for reading the file content
    FileInputStream fis = new FileInputStream(file);

    //Create byte array to read data in chunks
    byte[] byteArray = new byte[1024];
    int bytesCount = 0;

    //Read file data and update in message digest
    while ((bytesCount = fis.read(byteArray)) != -1) {
        digest.update(byteArray, 0, bytesCount);
    };

    //close the stream; We don't need it now.
    fis.close();

    //Get the hash's bytes
    byte[] bytes = digest.digest();

    //This bytes[] has bytes in decimal format;
    //Convert it to hexadecimal format
    StringBuilder sb = new StringBuilder();
    for(int i=0; i< bytes.length ;i++)
    {
        sb.append(Integer.toString((bytes[i] & 0xff) + 0x100, 16).substring(1));
    }

    //return complete hash
   return sb.toString();
}

/**
 * Test if FASTQC runs with single-end data
 */
workflow test_single_end {
    input_files = Channel.fromPath("${baseDir}/input/test_single_end.fastq.gz")
                    .map {f -> [f.name.replace(".fastq.gz", ""), true, f]}
    FASTQC(input_files)

    // test that the output looks as expected
    FASTQC.out.html.map { name, is_single_end, html_file ->
        html_hash = getFileChecksum(MessageDigest.getInstance("MD5"), new File("${html_file}"));

        assert name == "test_single_end"
        assert is_single_end == true
        assert html_file.getName() == "test_single_end_fastqc.html"
        // Hash seems to vary between local runs and GitHub Actions
        // TODO: Might be solved when using Docker for tests?
        // assert html_hash == "8ed68442ebb5b9706bf79b4f66701e15"
    }
    FASTQC.out.zip.map { name, is_single_end, zip_file ->
        // NOTE: output zip files do not have a consistent hash
        assert name == "test_single_end"
        assert is_single_end == true
        assert zip_file.getName() == "test_single_end_fastqc.zip"
    }
}

/**
 * Test if FASTQC runs with paired end data
 */
workflow test_paired_end {
    input_files = Channel.fromFilePairs("input/test_R{1,2}.fastq.gz")
                    .map {f -> [f[0], false, f[1]]}
    FASTQC(input_files)

    // test that the output looks as expected
    FASTQC.out.html.map { name, is_single_end, html_files ->
        html_r1 = html_files[0]
        html_r2 = html_files[1]

        html_r1_hash = getFileChecksum(MessageDigest.getInstance("MD5"), new File("${html_r1}"));
        html_r2_hash = getFileChecksum(MessageDigest.getInstance("MD5"), new File("${html_r2}"));

        assert name == "test_R"
        assert is_single_end == false
        assert html_r1.getName() == "test_R_1_fastqc.html"
        assert html_r2.getName() == "test_R_2_fastqc.html"
        assert html_r1_hash == "082c13ce7163ea0f52a66b83cb57b0f0"
        assert html_r2_hash == "4ff04ec8da77e3af512f03b8c09a9e04"
    }
    FASTQC.out.zip.map { name, is_single_end, zip_files ->
        zip_r1 = zip_files[0]
        zip_r2 = zip_files[1]
        // NOTE: output zip files do not have a consistent hash

        assert name == "test_R"
        assert is_single_end == false
        assert zip_r1.getName() == "test_R_1_fastqc.zip"
        assert zip_r2.getName() == "test_R_2_fastqc.zip"
    }
}

workflow {
    test_single_end()
    test_paired_end()
}
