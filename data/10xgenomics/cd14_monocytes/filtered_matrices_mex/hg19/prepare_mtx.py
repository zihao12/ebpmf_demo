## Input: Y.csv, Yhat.csv
## output: train.mtx, test.mtx, validation.mtx

from scipy import sparse, io
import numpy as np
import os
np.random.seed(123)

Y = np.genfromtxt("Y.csv",delimiter = ",")
Yhat = np.genfromtxt("Yhat.csv",delimiter = ",")

n_val = min(500, int(0.1*Y.shape[0]))
idx = np.random.randint(Y.shape[0], size = n_val)
Yval = Yhat[idx,:]

io.mmwrite("train.mtx", sparse.csc_matrix(Y.astype(int)))
io.mmwrite("test.mtx", sparse.csc_matrix(Yhat.astype(int)))
io.mmwrite("validation.mtx", sparse.csc_matrix(Yval.astype(int)))

os.system("tail -n +3 train.mtx > train.tsv")
os.system("tail -n +3 test.mtx > test.tsv")
os.system("tail -n +3 validation.mtx > validation.tsv")

