For development, ensure you are running a docker instance hosting a postgressql DB which contains previously reported pipestat results:

Docker command:
```
sudo docker run --rm -it -e POSTGRES_USER=postgres -e POSTGRES_PASSWORD=pipestat-password -e POSTGRES_DB=pipestat-test -p 5432:5432 postgres
```

To run with reload on:
`cd pipestat/pipestatreader`
`export PIPESTAT_CONFIG="path/to/pipestatconfigfile.yaml"`
`uvicorn reader:app --reload --port 8000`

To run via pipestat cli:
`pipestat serve --config "path/to/pipestatconfigfile.yaml" --host "0.0.0.0" --port 8000`
