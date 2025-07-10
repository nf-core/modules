FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:3da7c4198dc77b50aeaf76f262ed0ac2ed324f87fba4b5b0f0bc0b4fbbf2ad93 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:f90174159add063f40041160baf4ead693302a285c12fa74e60154f1d7acbc16

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
