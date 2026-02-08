"""
RiboStructMapper - Map ribosome density data onto protein structures

A bioinformatics pipeline that integrates ribosome profiling data with protein 
structures to generate visualization-ready PDB files with density scores encoded 
in the B-factor column.
"""

__version__ = "1.0.0"
__author__ = "RiboStructMapper Team"

# Main package exports will be available to users
from ribostruct.core import parser, alignment, processor, injector

__all__ = ['parser', 'alignment', 'processor', 'injector']
