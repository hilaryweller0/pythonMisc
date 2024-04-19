a=2
assert a == 2, "hello"
assert a != 2, "a is not equal to 2"
assert expr, msg


try:
   func(para, meter)
   raise Exception
except exception:
   pass
   
import numpy as np

try:
    np.sqrt('a')
except TypeError:
   print 'Cannot take the square root of a string'

assert np.sqrt(1) == 1, "square root of 1 should be 1"
