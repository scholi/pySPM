import os
import sys

cmd = None

if len(sys.argv)>1:    
    cmd = sys.argv[1]
  
commands = {}

def script(func):
    __requires__ = 'pySPM'
    __import__('pkg_resources').run_script('pySPM', func)
    
cmds =  {
    'stability': dict(doc="Run a stability plotter tool in a GUI", script='stability.py'),
    'plotter':   dict(doc="Display a live graphical display of all the logged parameters (usually the emission current)", alias=['emission_current_plotter','parameters_plotter'], script = 'emission_current_plotter.py'),
    'timer': dict(doc="Display informations about the current measurement time. In particular this function estimate the remaining time from the measurement proportion and the elapsed time.", alias=['tof_timer','measurement_timer'], script='tof_timer.py')

for c in cmds:
    doc = cmds[c].get('doc', '')
    script = cmds[c].get('script', None)
    commands[c] = (script, doc)
    if cmds.has_key('alias'):
        for a in cmds[c]['alias']:
            commands[a] = (script, doc)
        
if cmd is None or cmd in ['help','--help','-h','?','-?']:
    print("Give as argument the command that pySPM should run. Valid commands are:")
    for c in cmds:
        print("{}: {} (aliases: {})".format(c, cmds[c].get('doc',''), ",".join(cmds[c].get('alias',[])))
elif cmd in commands:
    print("Running command \"{}\"".format(cmd))
    commands[cmd][0]()
else:
    print("Command \"{}\" not understood".format(cmd))
    print("Valid commands are:")
    for c in commands:
        print("{}: {}".format(c,commands[c][1]))