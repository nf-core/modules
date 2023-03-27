Except the sentieon bwaindex module, these nextflow sentieon modules need the sentieon license string (probably some url) to be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. The encoding can be done like this:

```
SENTIEON_LICENSE_BASE64=$(echo -n <sentieon_license_string> | base64 -w 0)
```

The nextflow secret is then set like this:

```
nextflow secret set SENTIEON_LICENSE_BASE64 $SENTIEON_LICENSE_BASE64
```
