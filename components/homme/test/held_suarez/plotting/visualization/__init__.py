# Append the root directory to the PYTHONPATH for future imports
import os, sys
root_dir: str = os.path.normpath( \
    os.path.join(os.path.dirname(__file__), os.pardir))
if root_dir not in sys.path:
    sys.path.append(root_dir)
    
# Standard Library Imports

# Third-Party Imports

# Local Library Imports

# Import from Local Source Files
from .plot_climatology import plot_climatology