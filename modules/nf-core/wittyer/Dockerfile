# syntax=docker/dockerfile:1-labs
FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:ee0d20fbc3dfc60ae2d9d9115a9b750c18576a47f24073119dd09d93e072cb89 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:e64ff42d964bcb32ccd4644ea5e4b114ffaf10d91df55c11ce895d670dce23d3

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
