FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:a574e6212c98ccf9cf2c6588201b1c71d93094f6f26c21375f75aa909a9e88f0 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:9916fb650c83b13231e4ec91d10e76a520209a0caeeeb46c507a9d98db8612a2

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
