import os
import datetime
import socket


a = 1
a.test = 1
print(a.__dir__())


def foo():
    return 1

foo.test = 1
print(dir(foo))
print(foo.__dict__)


