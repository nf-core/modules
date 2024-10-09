FROM mcr.microsoft.com/dotnet/sdk:6.0@sha256:ba7877fb3a171ad4326caf97d8b56304fe73e6e05da0992775804805f6f75451 as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:f6cf89e960d5f478b36ec19146dcdcd754c2369dbd0f9118c8fe444fff3dcf9d
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
