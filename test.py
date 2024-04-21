import os
import datetime
import socket


import time
import numpy as np
def calc_fib(n):
    if n <= 1:
        return n
    else:
        return calc_fib(n - 1) + calc_fib(n - 2)

var = 40
start = time.time()
print(f"{var}th number of Fibonacci is {calc_fib(var)}")
end = time.time()
print(f"time: {end - start} sec.")



