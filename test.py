from time       import *
from random     import *
from subprocess import *
import os

def randomlog(filename):
    f = open(filename, 'w')
    for _ in range(10):
        sleep(1)
        f.write("%0d\n" % randint(-1000,1000))
    f.close


pipename = "./pipe1"

os.mkfifo(pipename)


if(os.fork() == 0):
    randomlog(pipename)
else:
    p = open(pipename)
    print(p.readlines())
    p.close();
    os.remove(pipename)