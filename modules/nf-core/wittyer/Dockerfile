FROM mcr.microsoft.com/dotnet/sdk:6.0@sha256:10655413946012e357306f0c14fd53341d6cbe391072b0d9a3c09c72eb6774ee as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:e21b5579c1c30e8245039f2d72e8b4d0a80f617cc49d3700e9e9992e07dde2bf
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
