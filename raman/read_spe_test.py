import numpy as np
import pandas as pd 
from tkinter import filedialog
import matplotlib.pyplot as plt
import os 
import glob
import csv
import itertools
from pyspec.ccd.files import PrincetonSPEFile

initialdir='/mnt/c/Users/roflc/Downloads/Solidstate/Raman/514.5nm/Redux/BLANK4R.ASC' #assign a directory 

f=PrincetonSPEFile(initialdir)
f.getData()