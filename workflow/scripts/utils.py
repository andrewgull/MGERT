import pandas as pd

def read_genome_paths(file_path: str) -> pd.DataFrame:
    """
    read paths to genome assembly files from a TSV file.
    arguments:
        file_path (str): Path to the TSV file containing genome assembly paths.
    returns:
        list of str: A list of genome assembly file paths.
    """
    try:
        df = pd.read_csv(file_path, sep="\t", header=None, names=["name", "path"])
        return df["path"].tolist()
    except Exception as e:
        raise ValueError(f"Error reading assembly paths from {file_path}: {e}")

def read_genome_names(file_path: str) -> pd.DataFrame:
    """
    read genome names from a TSV file.
    arguments:
        file_path (str): Path to the TSV file containing genome names.
    returns:
        list of str: A list of genome names.
    """
    try:
        df = pd.read_csv(file_path, sep="\t", header=None, names=["name"])
        return df["name"].tolist()
    except Exception as e:
        raise ValueError(f"Error reading genome names from {file_path}: {e}")

def get_path_using_name(name: str, file_path: str) -> str:
    """
    Get the path of a genome assembly using its name.
    arguments:
        name (str): The name of the genome.
        file_path (str): Path to the TSV file containing genome names and paths.
    returns:
        str: The path of the genome assembly.
    """
    df = pd.read_csv(file_path, sep="\t", header=None, names=["name", "path"])
    row = df[df["name"] == name]
    if not row.empty:
        return row["path"].values[0]
    else:
        raise ValueError(f"Genome name '{name}' not found in {file_path}.")
