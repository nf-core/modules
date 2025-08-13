FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:ae000be75dac94fc40e00f0eee903289e985995cc06dac3937469254ce5b60b6 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:52011e65ae82d61566a8a028ef72e0a77dbbab6c9233212fc28ad67bc55f1354

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
