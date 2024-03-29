import sys

cmd = None

if len(sys.argv) > 1:
    cmd = sys.argv[1]

commands = {}


def run_script(func):
    import importlib

    mod = importlib.import_module("pySPM.tools." + func)
    sys.argv = [sys.argv[0], *list(sys.argv[2:])]
    mod.main()


cmds = {
    "stability": {
        "doc": "Run a stability plotter tool in a GUI",
        "script": "stability",
    },
    "plotter": {
        "doc": "Display a live graphical display of all the logged parameters (usually the emission current)",
        "alias": ["emission_current_plotter", "parameters_plotter"],
        "script": "emission_current_plotter",
    },
    "timer": {
        "doc": "Display informations about the current measurement time. In particular this function estimate the remaining time from the measurement proportion and the elapsed time.",
        "alias": ["tof_timer", "measurement_timer"],
        "script": "tof_timer",
    },
    "spectra": {"doc": "Spectra plotter and mass calibration GUI", "script": "spectra"},
}

for c in cmds:
    doc = cmds[c].get("doc", "")
    script = cmds[c].get("script", None)
    commands[c] = (script, doc)
    if "alias" in cmds[c]:
        for a in cmds[c]["alias"]:
            commands[a] = (script, doc)

if cmd is None or cmd in ["help", "--help", "-h", "?", "-?"]:
    print("Give as argument the command that pySPM should run. Valid commands are:")
    for c in cmds:
        print(
            "{}: {} (aliases: {})".format(
                c, cmds[c].get("doc", ""), ",".join(cmds[c].get("alias", []))
            )
        )
elif cmd in commands:
    print(f'Running command "{cmd}"')
    run_script(commands[cmd][0])
else:
    print(f'Command "{cmd}" not understood')
    print("Valid commands are:")
    for c in commands:
        print(f"{c}: {commands[c][1]}")
