FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:51060ced0e5dc1282a587a092783132ba8b0ec25150d3bc61e99bf0ef2d62593 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:c8029dff6ec67e25907cb5fa4c7e26605cc1784606d8146c920a782b64dfa9e6

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
