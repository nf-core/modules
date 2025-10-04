FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:123b43e4d9775451d8ed63af324a42132707b4edb14770bcdf5c85cf55bc45f1 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:6e6ccf450a40a99748e93bf18dd6739893881f5dd831ae23b586565749b1643d

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
