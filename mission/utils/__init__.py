# This file is created on 2024/10/16 for the purpose of reproducibility.


import pathlib

def get_project_root() -> pathlib.Path:
    """
    Get the path to the root of the project

    :return: The path to the root of the project
    """
    return pathlib.Path(__file__).parent.parent.parent

def get_data_dir() -> pathlib.Path:
    """
    Get the path to the data directory

    :return: The path to the data directory
    """
    project_root = get_project_root()
    return project_root / "data"

def get_species_genome_dir(species: str) -> pathlib.Path:
    """
    Get the path to the genome directory for a given species

    :param species: The species to get the genome directory for; one of "human", "chicken", "fruitfly", or "worm"
    :return: The path to the genome directory for the given species
    """
    data_dir = get_data_dir()
    if species == "human":
        return data_dir / "hg38" / "chroms"
    elif species == "chicken":
        return data_dir / "galGal3"
    elif species == "fruitfly":
        return data_dir / "dm3"
    else:
        return data_dir / "ce10"
