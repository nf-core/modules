# How to use the sentieon modules

Except the sentieon bwaindex module, these nextflow sentieon modules require a license. Sentieon supply license in the form of a string-value (a url) or a file. It should be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. If a license string (url) is supplied, then the nextflow secret should be set like this:

```bash
nextflow secret set SENTIEON_LICENSE_BASE64 $(echo -n <sentieon_license_string> | base64 -w 0)
```

If a license file is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
```
