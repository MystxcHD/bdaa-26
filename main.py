import pandas as pd

file_path = "data/GSE33000_series_matrix (1).txt"

data = pd.read_csv(
    file_path,
    sep="\t",
    skiprows=88,
    nrows=39369 - 88,
    index_col=0
)

print(data.head())
print(data.shape)
data = data.T

print(data.head())
print(data.shape)