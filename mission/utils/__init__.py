# This file is created on 2024/10/16 for the purpose of reproducibility.


import pathlib

def get_project_root():
    return pathlib.Path(__file__).parent.parent.parent

def get_data_dir():
    project_root = get_project_root()
    return project_root / "data"
