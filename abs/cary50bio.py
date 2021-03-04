from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

directory = '/mnt/c/users/roflc/Desktop/Python Test Data/'
filename = 'Nickel(II)Sulfate.csv'
filepath = directory + filename

data=pd.read_csv(filepath)
print(data.head(5))