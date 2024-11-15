FROM mcr.microsoft.com/dotnet/sdk:6.0@sha256:aca95d87aec3d66c4874d25a06845d1958de588de747c248f735f05db0c396f3 as builder

WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:93be1090f1a66815a580b844aedd33b5f85e4821209c606217de95d73a10d6d6

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
