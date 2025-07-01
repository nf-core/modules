FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:b768b444028d3c531de90a356836047e48658cd1e26ba07a539a6f1a052a35d9 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:267364a9e6d272b64dcb490ecdddbebaee87b8ed13e911e2e2cb7d14d51e085a

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
