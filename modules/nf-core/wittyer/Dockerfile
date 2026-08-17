# syntax=docker/dockerfile:1-labs
FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:35048e3a81e6a07c316e7bbbd80d80d2ba705fe5f23a8ed42b6638c8f4c20d30 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:ee9e6309cef467e134056f9115b31fe3a43ef2959e5b1f42bce0cc97f1b3db3f

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
