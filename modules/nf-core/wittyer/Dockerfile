# syntax=docker/dockerfile:1-labs
FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:6416205d95f47e66d4292310710f7d738a4ad8f499c8d2370cf22d29a3c7d4b3 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:949545c232b0e78c15177f6b60f64d7d6cd835a98279fac506d9a4512a5e0151

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
