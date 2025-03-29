# Updating the docker container and making a new module release

Space Ranger is a commercial tool by 10X Genomics. The container provided for the spaceranger nf-core module is not provided nor supported by 10x Genomics. Updating the Space Ranger version in the container and pushing the update to Dockerhub needs to be done manually.

1. Navigate to the [Space Ranger download page](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest) and get a link for the Space Ranger downloads

2. Edit the Dockerfile: update the Space Ranger version in this line:

```diff
- ARG SPACERANGER_VER="3.1.2"
+ ARG SPACERANGER_VER="<New version>"
- ARG SPACERANGER_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.2.tar.gz?Expires=1732608367&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=S4jcGCs5H0wLFfREhSc~AfuPIOXE~PW96kX9o2EyxZnmn1goUecgfEWTj67fL1DnZJBIa06kFdUdNpsMn2ustiSWOfXEgjlYQnxIcNnQIiEBGoQTHWphvx3LEQ6wtZnkWS80P6IcE0HJkIsgy04t6Sohih5cxY4jgytYsrAfZDYr5G3KKFwTfCKmhzMaXqW635yPbyQ8xEcQHK0QwviAx8-EFq-PE8UzC4QgUKi2MW-ivcfZkSDSfF8C3s7SgwDXIGIWv52mzeszenxMjN4KrWQotZ7ZpktzI0Vfpz0dNC17dQeDQUHj4LuNYbdh3RqsPKtqu3wjCe2Q7KiyoWnmaw__" \"
+ ARG SPACERANGER_URL="<Your-unique-download-link-that-expires"
- ARG SPACERANGER_SHA256="2566b24f29829b39f3add112a674990b1c54ae2fbe7ccb50a4c7dce9ccf152e6"
+ ARG SPACERANGER_SHA256="<Manual SHA256 of tar.gz>"
```

3. Push the changes and the Dockerfile should be built and uploaded to `quay.io/nf-core/modules/spaceranger`!
