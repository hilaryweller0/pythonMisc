#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
import pandas as pd
import datetime as dt
import math

# Test math
assert math.factorial(3) == 6

# Test datetime
date = dt.datetime.strptime('20200928', '%Y%m%d')
assert date.weekday() == 0

# Test numpy
a = np.linspace(0,1,3)
assert a[1] == 0.5

# Test linanlg
A = np.zeros([2,2])
A[0,0] = 2
A[1,1] = 2
b = la.solve(A, np.ones(2))
assert b[0] == 0.5

# Test pandas write and read
df1 = pd.DataFrame([['a', 'b'], ['c', 'd']],
                   index=['row 1', 'row 2'],
                   columns=['col 1', 'col 2'])

df1.to_excel("output.xlsx")  

df2 = pd.read_excel('output.xlsx', index_col=0)
assert df2.values[0,0] == 'a'

# Test matplotlib
x = np.linspace(0,2*np.pi, 41)
y = np.sin(x)
plt.plot(x,y)
plt.show()

