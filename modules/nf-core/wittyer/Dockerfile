FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:584cdb686f8ff7d2ffd636581b9bba8c9de38df1e751e9275438cddecdfbbcd2 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:54fb43f9f282a6725e0b6462fe7f5152d9611404af8cffb4613a31f714551078

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
