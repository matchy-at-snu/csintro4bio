import os
import sys

# Get the absolute path of the current script
current_script_path = os.path.abspath(__file__)

# Navigate up to the parent directory of 'mission'
project_root = os.path.dirname(os.path.dirname(current_script_path))

# Add the project root to the Python path
sys.path.insert(0, project_root)

# Now you can import from mission.utils
from mission.utils import *

# Your mission code goes here
