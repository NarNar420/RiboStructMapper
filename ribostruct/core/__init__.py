"""
Core processing modules for RiboStructMapper.

This package contains the fundamental modules for:
- Parsing genomic data, PDB structures, and density files
- Sequence translation and alignment
- Density aggregation and offset processing
- B-factor injection into PDB structures
"""

from ribostruct.core import parser, alignment, processor, injector

__all__ = ['parser', 'alignment', 'processor', 'injector']
