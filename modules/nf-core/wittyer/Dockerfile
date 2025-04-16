FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:7d24e90a392e88eb56093e4eb325ff883ad609382a55d42f17fd557b997022ca as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:1e5eb0ed94ca96a34a914456db80e48bd1bb7bc3e3c8eda5e2c3d89c153c3081

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
