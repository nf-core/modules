FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:90872f8e7f1fd2b93989b81fb7f152c3bef4fe817470a3227abaa18c873dba60 as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:14529413f5bcd9f11d3fb230abe5bcc6118ca27cc3a3d3ee042dc36081a9569f

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
