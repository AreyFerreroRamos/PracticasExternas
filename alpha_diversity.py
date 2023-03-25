# Author: Arey Ferrero Ramos.
# Date: March 24, 2023. Version: 1
# Description: This script calculates the alpha diversity of a species of vertebrate from those used in the study.
#   Parameters:
#       -A file that contains a table with the number of species per genus in every individual of that species of
#        vertebrate used in that study.
#   Output:
#       -The alpha diversity of the species of vertebrate.

import pandas as pd
import sys
import os

if len(sys.argv) == 2:
    if os.path.isfile(sys.argv[1]):
        file = pd.read_csv(sys.argv[1], sep=' ')
        print(file)
    else:
        print("Error: The input must be a file corresponding to a species of vertebrate.")
else:
    print("Error: The number of parameters is incorrect. The only parameter has to be a file.")