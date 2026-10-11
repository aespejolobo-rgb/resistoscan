"""resistoscan: alignment-free screening of antimicrobial resistance genes.

Detects and quantifies antimicrobial resistance genes (ARGs) in clinical
isolate assemblies and environmental shotgun metagenomes, supporting One
Health surveillance of antimicrobial resistance.
"""
__version__ = "0.1.0"

from .index import KmerIndex  # noqa: E402
from .screen import ReadScreenParams, screen_assembly, screen_reads  # noqa: E402

__all__ = ["KmerIndex", "ReadScreenParams", "screen_assembly", "screen_reads", "__version__"]
