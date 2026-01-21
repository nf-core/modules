#!/usr/bin/env python3

import platform

import metaspace
import pandas as pd
from metaspace import SMInstance


def download_dataset_results(dataset_id, database=None, version=None):
    """
    Download the results of a specified dataset and save them as a CSV file.

    Args:
        dataset_id (str): The ID of the dataset to process.
            The ID is the last part of the dataset URL.
            For example, in the URL `https://metaspace2020.org/dataset/2022-08-05_17h28m56s`,
            the `dataset_id` is `2022-08-05_17h28m56s`.
        database (str, optional): The name of the database to filter. If not provided, all databases will be included.
        version (str, optional): The version of the database to filter. If not provided, all versions will be included.

    Returns:
        pd.DataFrame: A DataFrame containing the results.

    Behavior:
        - If only `database` is provided and it has multiple versions, the results will be combined into one file named `{database}_all_versions.csv`.
        - If the `database` has only one version, the file will be named `{database}_{version}.csv`.
        - If no `database` is provided, it will download all the annotations file from all database and the file will be named `all_databases.csv`.
    """
    # Retrieve the dataset object
    print(f"Processing dataset_id: {dataset_id}, database: {database}, version: {version}")
    # print("Initializing SMInstance...")
    sm = SMInstance()
    print(f"Fetching dataset: {dataset_id}")
    try:
        dataset = sm.dataset(id=dataset_id)
    except Exception:
        print(f"❌ {dataset_id} Dataset not found or inaccessible. Skipping this dataset.")
        return None

    # Step 1: Generate annotation statistics table
    # print("Fetching database details...")
    databases = dataset.database_details
    # Extract database names and versions
    database_info = []
    for db in databases:
        db_str = str(db)
        name_version = db_str.split(":")[1:3]
        name_version = [nv.strip(">") for nv in name_version]  # Remove trailing '>'
        database_info.append(" ".join(name_version))

    print("Available databases and versions:")
    for db in database_info:
        print(f"  - {db}")

    # Check if the multiple versions database exists
    has_multiple_versions = False
    # Check if the specified database exists
    if database:
        matched_databases = [db for db in database_info if db.startswith(database)]
        if not matched_databases:
            print(f"❌ {dataset_id} could not find database: {database}")
            return None
        elif len(matched_databases) > 1 and version is None:
            has_multiple_versions = True
            print(f"Warning: Multiple versions found for database {database}.")
            print("Available versions:")
            for db in matched_databases:
                print(f"  - {db}")

    # Initialize the results DataFrame
    results_df = pd.DataFrame()

    # Process the specified database names and versions
    for db in database_info:
        name, ver = db.split()
        if (database is None or database == name) and (version is None or version == ver):
            print(f"Processing database: {name}, version: {ver}")
            result = dataset.results(database=(name, ver))
            result_df = pd.DataFrame(result)
            result_df["Database"] = name
            result_df["Version"] = ver
            results_df = pd.concat([results_df, result_df], ignore_index=True)

    if results_df.empty:
        print(f"❌ {dataset_id} has no annotation data in database: {database}.")
        return None

    # Extract the content of '+adduct' and '-adduct''
    if "ion" in results_df.columns:
        results_df["Adduct"] = results_df["ion"].str.extract(r"([+-][A-Za-z0-9]+)")
    # Rename the `intensity` column to `maxIntensity`
    if "intensity" in results_df.columns:
        results_df.rename(columns={"intensity": "maxIntensity"}, inplace=True)

    # Dynamically generate the CSV file name
    if database:
        matched_versions = [db.split()[1] for db in database_info if db.startswith(database)]
        if len(matched_versions) > 1:
            filename = f"{dataset_id}_{database}_all_versions.csv"
        else:
            filename = f"{dataset_id}_{database}_{matched_versions[0]}.csv"
    else:
        filename = f"{dataset_id}_all_databases.csv"

    # Save the results to a CSV file
    results_df.to_csv(filename, index=False)
    if has_multiple_versions:
        print(f"⚠️ {dataset_id} has multiple {database} version. All the version saved to {filename}")
    elif database:
        print(f"✅ {dataset_id} with {database} database are saved to {filename}")
    else:
        print(f"✅ {dataset_id} with all database are saved to {filename}")

    return results_df


dataset_id = "${dataset_id}"

if dataset_id == "null" or dataset_id is None:
    raise ValueError("Error: input csv or datasets contains an entry with missing dataset_id.")

database = "${database}" if "${database}" != "null" else None
version = "${version}" if "${version}" != "null" else None

result = download_dataset_results(dataset_id, database, version)


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "    " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "metaspace": metaspace.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
print("versions.yml file has been generated.")
