import os
import sys

current_script_path = os.path.abspath(__file__)

project_root = os.path.dirname(os.path.dirname(current_script_path))

sys.path.insert(0, project_root)

from mission.utils import *
