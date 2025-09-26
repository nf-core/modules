FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:bb42ae2c058609d1746baf24fe6864ecab0686711dfca1f4b7a99e367ab17162 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:5040af5025cf2e12b1c4eb04c033342c39fe654bbd9395d24cc2635ec196c7f4

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
