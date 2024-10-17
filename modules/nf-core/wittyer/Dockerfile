FROM mcr.microsoft.com/dotnet/sdk:6.0@sha256:3e7491422fb3a14f0bfe98fcb58afc1d477ee6ee4aac8dff8557332140ac9cf1 as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:737d5276dcae2a4ba8cf5b860f0b4147f5f1fb6749e78232c2aaba2e825e4061
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
