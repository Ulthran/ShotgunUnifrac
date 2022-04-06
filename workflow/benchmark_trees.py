import unifrac
import numpy as np
from skbio import DistanceMatrix

dms: list = []
for fp in snakemake.input:
    dms.append(unifrac.unweighted_fp32(phylogeny=fp))

### Considerations
#
# Distance matrices might not be same size
# Rows/columns might be unordered
# 