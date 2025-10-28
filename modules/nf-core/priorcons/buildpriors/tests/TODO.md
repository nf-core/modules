# TODO: Update test file

## Tasks to complete for `main.nf.test`

1. **Update test inputs** to match the module's actual inputs:
   - `input[0]`: tuple with meta and alignment FASTA file
   - `input[1]`: reference sequence ID (string)
   - Currently uses BAM file - need to change to FASTA alignment

2. **Update test data paths**:
   - Find appropriate test alignment FASTA in `params.modules_testdata_base_path`
   - Or add new test data if alignment FASTA doesn't exist

3. **Update test names**:
   - Change from "sarscov2 - bam" to something like "sarscov2 - fasta alignment"
   - Keep the pattern: `test_name`, `test_name - stub`

4. **Optional: Add ext.args example**:

   ```groovy
   when {
       process {
           """
           input[0] = [
               [ id:'test' ],
               file('path/to/alignment.fasta')
           ]
           input[1] = 'MN908947.3'  // reference ID

           // Optional: demonstrate custom parameters
           ext.args = '--win 50 --overlap 25'
           """
       }
   }
   ```

5. **Update output assertions**:
   - Verify output is `.parquet` file (not `.bam`)
   - Update snapshot expectations

6. **Run tests** after updates:
   ```bash
   nf-core modules test priorcons/buildpriors
   ```

## Reference

- nf-test docs: https://nf-co.re/docs/contributing/modules
- Chained modules: https://nf-co.re/docs/contributing/modules#steps-for-creating-nf-test-for-chained-modules
