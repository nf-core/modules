FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:ecdbce3ac370d065b4e34bdb367ed8ac10e24d2af630716f9828bf24b360691c as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:ba038432e8b654fa13461de3b3a03089a23c70e05ea0aac29f6389a4b654e768
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
