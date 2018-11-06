from .win32_helper import *
import psutil

class Fpanel:
    def __init__(self):
        self.PID = getPID('fpanel5.exe')
        if self.PID is not None:
            self.proc = psutil.Process(self.PID)
            self.hwnd = findWindow("Fpanel")[0]
        else:
            raise Exception("Fpanel is not running")
        
    def isInstrumentPanelOpen(self):
        try:
            hwnd = findTopWindow(wantedClass="InstrumentPanel.36575160")
            return hwnd
        except:
            return False

    def getLogFile(self):
        logs = []
        for f in self.proc.open_files():
            if f.path.startswith("C:\\Program") or f.path.startswith("C:\\Users") or f.path.startswith("C:\\Windows"):
                continue
            if f.path.endswith(".txt") or f.path.endswith(".log"):
                logs.append(f.path)
        if len(logs)==0:
            log = None
        elif len(logs)==1:
            log = logs[0]
        else:
            logs2 = [x for x in logs if 'log' in x]
            if len(logs)==1:
                log = logs2[0]
            else:
                log = logs[0]
        return log

    def openInstrumentPanel(self):
        controls = []
        for x in findControls(self.hwnd, wantedClass="Button"):
            bbox = getBBox(x)
            w = bbox[2]-bbox[0]
            h = bbox[3]-bbox[1]
            if w!=28 or h!=28:
                continue
            controls.append([x]+bbox)
        Y = sorted(list(set([x[2] for x in controls])))
        hInstrumentButton = [x[0] for x in controls if x[2]==Y[3]][0]
        clickButton(hInstrumentButton)
        clickButton(hInstrumentButton)

    def openTab(self, i):
        if not self.isInstrumentPanelOpen():
            self.openInstrumentPanel()
        hwnd = findTopWindow(wantedClass="InstrumentPanel.36575160")
        hwnd = selectVert(selectHoriz(findControls(hwnd, wantedClass="Button"), 0), i)[0]
        clickButton(hwnd)
        clickButton(hwnd)
        
    def openLogWindow(self):
        if not self.isLogWindowOpen():
            self.openTab(8)
            hwnd = findTopWindow(wantedClass="InstrumentPanel.36575160")
            hBt = findControl(hwnd, wantedText="Log...", wantedClass="Button")
            clickButton(hBt)
            clickButton(hBt)
        return findTopWindow(wantedClass="LogFrame.259931064")
        
    def isLogWindowOpen(self):
        hwnd = findWindow("Log Settings")
        if len(hwnd)==0:
            return False
        return hwnd[0]

    def populateLogFields(self):
        hlog = self.openLogWindow()
        hwnd = findTopWindow(wantedClass="InstrumentPanel.36575160")
        hwnds = findControls(hwnd, wantedClass="PythonCustomControl.Meter.36566376")
        if True:
            hAdd = findControl(hlog, "Add", "Button")
            clickButton(hAdd)
            self.openTab(1)
            hEC = selectVert(hwnds, 1)[0]
            bbox = getBBox(hEC)
            x = (bbox[0]+bbox[2])//2
            y = (bbox[1]+bbox[3])//2
            click(x, y)