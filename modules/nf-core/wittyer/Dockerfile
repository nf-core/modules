FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:2212719b38aa69ec476e181e41c4a4bff05cb3dd4904e0a999c7bb50a614e318 as builder
WORKDIR /src
COPY . /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:566d713e9f697ed7a1ff54fc467642ce5764bd06399cb204ecae478d485883f8
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
