#!/usr/bin/env python
#
#  multiprocessing_test1.py
#  
#
#  Created by Danny Jacobs on 3/23/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,multiprocessing as mp,os
from os import environ
from multiprocessing import Process,Lock,Array,Value,Pool
import os
import time

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('--ex',type='int',default=1,
   help="Example number to do.")
opts, args = o.parse_args(sys.argv[1:])


if opts.ex==0:
    def my_fork():
        child_pid = os.fork()
        if child_pid == 0:
            print "Child Process: PID# %s" % os.getpid()
        else:
            print "Parent Process: PID# %s" % os.getpid()

    if __name__ == "__main__":
        my_fork()
elif opts.ex==1:
    def my_fork():
        environ['FOO']="baz"
        print "FOO environmental variable set to: %s" % environ['FOO']
        environ['FOO']="bar"
        print "FOO environmental variable changed to: %s" % environ['FOO']
        child_pid = os.fork()
        if child_pid == 0:
            print "Child Process: PID# %s" % os.getpid()
            print "Child FOO environmental variable == %s" % environ['FOO']
        else:
            print "Parent Process: PID# %s" % os.getpid()
            print "Parent FOO environmental variable == %s" % environ['FOO']

    if __name__ == "__main__":
        my_fork()
elif opts.ex==2:
    def sleeper(name, seconds):
       print 'starting child process with id: ', os.getpid()
       print 'parent process:', os.getppid()
       print 'sleeping for %s ' % seconds
       time.sleep(seconds)
       print "Done sleeping"


    if __name__ == '__main__':
       print "in parent process (id %s)" % os.getpid()
       p = Process(target=sleeper, args=('bob', 5))
       p.start()
       print "in parent process after child process start"
       print "parent process about to join child process"
       p.join()
       print "in parent process after child process join" 
       print "parent process exiting with id ", os.getpid()
       print "The parent's parent process:", os.getppid()
elif opts.ex==3:
    def printman(number,seconds):
        for i in range(number):
            print os.getpid(),
            time.sleep(seconds)
            print "Done!"
    if __name__=='__main__':
        p = Process(target=printman,args=(10,1))
        print "in parent process (id %s)" % os.getpid()        
        print "starting printman process"
        p.start()
        print "something happening after start()"
        print "sleeping for a second"
        time.sleep(1)
        print "joining up to printer process"
        p.join()
        print "main here: going to sleep"
        time.sleep(1)
        print "waking up!"
        print "EOL"
elif opts.ex==4:
    def f(l, i):
        l.acquire()
        print 'hello world', i
        l.release()

    if __name__ == '__main__':
        lock = Lock()

        for num in range(10):
            Process(target=f, args=(lock, num)).start()

elif opts.ex==5:
    print "testing shared memory"
    def f(n, a):
        n.value = 3.1415927
        for i in range(len(a)):
            a[i] = -a[i]

    if __name__ == '__main__':
        num = Value('d', 0.0)
        arr = Array('f', n.array(range(10)))

        p = Process(target=f, args=(num, arr))
        p.start()
        p.join()

    print num.value
    print arr[:]
elif opts.ex==6:
    print "Pools"
    def f(x):
        return x*x
    if __name__=='__main__':
        pool = Pool(processes=4)
        result = pool.apply_async(f,[10])
        print pool.map(f,range(10))
        print result.get(timeout=1)
       
