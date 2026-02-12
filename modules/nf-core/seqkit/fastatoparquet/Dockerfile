FROM debian

RUN apt-get update && apt-get install -y \
    wget unzip seqkit \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/duckdb/duckdb/releases/download/v1.0.0/duckdb_cli-linux-aarch64.zip \
    && unzip duckdb_cli-linux-aarch64.zip -d /usr/local/bin \
    && rm duckdb_cli-linux-aarch64.zip

RUN duckdb -c "INSTALL parquet;"; \
    echo "LOAD parquet;" >> $HOME/.duckdbrc;
