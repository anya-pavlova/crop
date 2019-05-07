import pandas as pd
import sys
import os

file_path = sys.argv[1]

df = pd.read_csv(file_path)

df_to_save = pd.DataFrame({"mzmed":df["mz"]}, index = df.index)
filename_without_extension = os.path.splitext(os.path.basename(file_path))[0]
df_to_save.to_csv(os.path.join("./input", filename_without_extension+".tsv"), sep="\t")