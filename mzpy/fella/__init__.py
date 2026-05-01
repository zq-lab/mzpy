"""
fella-py: Python port of the FELLA package for metabolomics enrichment.

rewrite from FELLA, a Bioconductor package for R, 
    by  Sergio Picart-Armada ,
          Francesc Fernandez-Albert ,
          Alexandre Perera-Lluna.
"""



from .fella_class import Fella
from .build_dataset import build_dataset
