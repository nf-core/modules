FROM mcr.microsoft.com/dotnet/sdk:9.0@sha256:9b0a4330cb3dac23ebd6df76ab4211ec5903907ad2c1ccde16a010bf25f8dfde as builder

ADD https://github.com/Illumina/witty.er.git#v0.5.2 /src
WORKDIR /src
RUN cd Wittyer \
    && dotnet publish -f net6.0 -r linux-x64 -c Release -o /output \
    && chmod +x /output/Wittyer

FROM mcr.microsoft.com/dotnet/runtime:9.0@sha256:17de18d43822c0379f9779ecaf5e8886b73b413b5009251a277c08ff7be47027

LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix libunwind8 openssl procps
COPY --from=builder /output /opt/Wittyer
