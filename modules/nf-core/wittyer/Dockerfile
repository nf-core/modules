FROM mcr.microsoft.com/dotnet/sdk:6.0@sha256:6df1177e48b55272316d08f19cb383483af82aca5cdc67a76c414bc200847624 as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:6.0@sha256:437cda84bdce26ce074d88b63abeec567c7226d73e8b911605077468e1d5c8d5
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
