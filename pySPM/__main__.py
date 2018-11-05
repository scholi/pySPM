import os
import sys

cmd = None

if len(sys.argv)>1:    
    cmd = sys.argv[1]
  
commands = {}

def script(func):
    
    global commands
    desc = func.__doc__.rstrip().lstrip()
    commands[func.__name__] = (func, desc)
    return func
    
@script
def stability():
    """
    Run a stability plotter tool in a GUI
    """
    __requires__ = 'pySPM>=0.2.9'
    print("Run stability plotter")
    __import__('pkg_resources').run_script('pySPM==0.2.9', 'stability.py')

if cmd is None or cmd in ['help','--help','-h','?','-?']:
    print("Give as argument the command that pySPM should run. Valid commands are:")
    for c in commands:
        print("{}: {}".format(c,commands[c][1]))
elif cmd in commands:
    print("Running command \"{}\"".format(cmd))
    commands[cmd][0]()
else:
    print("Command \"{}\" not understood".format(cmd))