FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:6c7713b7d107f743015b863f6bb028fca7430192140536779c0464f1ece899bc as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src

WORKDIR /src/Wittyer
RUN dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:de6d6aceb14f88c4ba6667f2fd87b0ce911800c5c6ed93b0906d21432b6a9386

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
