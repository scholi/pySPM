from ctypes import wintypes
import win32gui, win32api, win32con, ctypes
import struct
import array
import os, sys
import psutil
from win32gui import FindWindow, FindWindowEx, SendMessage, PyGetString
from win32con import WM_GETTEXT, WM_GETTEXTLENGTH

def getText(hwnd):
    buffer_len = SendMessage(hwnd, WM_GETTEXTLENGTH, 0, 0) + 1
    buffer = array.array('b', b'\x00\x00' * buffer_len)
    text_len = SendMessage(hwnd, WM_GETTEXT, buffer_len, buffer)
    text = PyGetString(buffer.buffer_info()[0], buffer_len - 1)
    return text

WNDENUMPROC = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)
user32 = ctypes.windll.user32
user32.EnumWindows.argtypes = [ WNDENUMPROC,wintypes.LPARAM]
user32.GetWindowTextLengthW.argtypes = [wintypes.HWND]
user32.GetWindowTextW.argtypes = [ wintypes.HWND, wintypes.LPWSTR,ctypes.c_int]

prototype = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, ctypes.POINTER(wintypes.RECT))
paramflags = (1, "hwnd"), (2, "lprect")
GetWindowRect = prototype(("GetWindowRect", user32), paramflags)

def clickButton(hwnd):
    _sendNotifyMessage(hwnd, win32con.BN_CLICKED)

def _getMultipleWindowValues(hwnd, getCountMessage, getValueMessage):
    result = []

    VALUE_LENGTH = 256

    valuecount = win32gui.SendMessage(hwnd, getCountMessage, 0, 0)
    for itemIndex in range(valuecount):
        valuebuffer = ctypes.create_unicode_buffer(VALUE_LENGTH)
        valueLength = win32gui.SendMessage(hwnd,
                                           getValueMessage,
                                           itemIndex,
                                           valuebuffer)
        result.append(valuebuffer[:valueLength])
    return result

def getListboxItems(hwnd):
    return _getMultipleWindowValues(hwnd,
                                     getCountMessage=win32con.LB_GETCOUNT,
                                     getValueMessage=win32con.LB_GETTEXT)

def findTopWindow(wantedText=None, wantedClass=None, selectionFunction=None):
    topWindows = findTopWindows(wantedText, wantedClass, selectionFunction)
    if topWindows:
        return topWindows[0]
    else:
        raise WinGuiAutoError("No top level window found for wantedText=" +
                               repr(wantedText) +
                               ", wantedClass=" +
                               repr(wantedClass) +
                               ", selectionFunction=" +
                               repr(selectionFunction))

def findTopWindows(wantedText=None, wantedClass=None, selectionFunction=None):
    results = []
    topWindows = []
    win32gui.EnumWindows(_windowEnumerationHandler, topWindows)
    for hwnd, windowText, windowClass in topWindows:
        if wantedText and not _normaliseText(wantedText) in _normaliseText(windowText):
            continue
        if wantedClass and not windowClass == wantedClass:
            continue
        if selectionFunction and not selectionFunction(hwnd):
            continue
        results.append(hwnd)
    return results

def findWindow(name='',C='',parent=0):
    R=[]
    def callbck(hwnd, lParam):
        title,c = getInfo(hwnd)
        if (name=='' or title.startswith(name)) and (C=='' or c==C):
            R.append(hwnd)
        return True
    user32.EnumChildWindows(parent,WNDENUMPROC(callbck),42)        
    return R

def listWin(hwnd, lParam):
    title,c = getInfo(hwnd)
    print(hwnd,title,c)
    return True

def getInfo(hwnd):
    length = user32.GetWindowTextLengthW(hwnd) + 1
    buffer = ctypes.create_unicode_buffer(length)
    user32.GetWindowTextW(hwnd, buffer, length)
    b2 = ctypes.create_unicode_buffer(200)
    user32.GetClassNameW(hwnd,b2,200);
    title=buffer.value
    c=b2.value
    return title,c

def showWin(parent,level=0):
    for x in findWindow(parent=parent):
        t,c = getInfo(x)
        print("\t"*level+hex(x)+" - "+t + " - " + c)
        showWin(x,level+1)

def getBBox(hwnd):
    b=GetWindowRect(hwnd)
    return [b.left,b.top,b.right,b.bottom]

def getPID(name):
    for proc in psutil.process_iter():
        if proc.name() == name:
            return proc.pid

def _windowEnumerationHandler(hwnd, resultList):
    '''Pass to win32gui.EnumWindows() to generate list of window handle,
    window text, window class tuples.'''
    resultList.append((hwnd,
                       win32gui.GetWindowText(hwnd),
                       win32gui.GetClassName(hwnd)))

def _normaliseText(controlText):
    '''Remove '&' characters, and lower case.
    Useful for matching control text.'''
    return controlText.lower().replace('&', '')

class WinGuiAutoError(Exception):
    pass

def findControl(topHwnd,
                wantedText=None,
                wantedClass=None,
                selectionFunction=None):
    controls = findControls(topHwnd,
                            wantedText=wantedText,
                            wantedClass=wantedClass,
                            selectionFunction=selectionFunction)
    if controls:
        return controls[0]
    else:
        raise WinGuiAutoError("No control found for topHwnd=" +
                               repr(topHwnd) +
                               ", wantedText=" +
                               repr(wantedText) +
                               ", wantedClass=" +
                               repr(wantedClass) +
                               ", selectionFunction=" +
                               repr(selectionFunction))

def findControls(topHwnd,
                 wantedText=None,
                 wantedClass=None,
                 selectionFunction=None):
    
    def searchChildWindows(currentHwnd):
        results = []
        childWindows = []
        try:
            win32gui.EnumChildWindows(currentHwnd,
                                      _windowEnumerationHandler,
                                      childWindows)
        except win32gui.error:
            # This seems to mean that the control *cannot* have child windows,
            # i.e. not a container.
            return
        for childHwnd, windowText, windowClass in childWindows:
            descendentMatchingHwnds = searchChildWindows(childHwnd)
            if descendentMatchingHwnds:
                results += descendentMatchingHwnds

            if wantedText and \
               not _normaliseText(wantedText) in _normaliseText(windowText):
                continue
            if wantedClass and \
               not windowClass == wantedClass:
                continue
            if selectionFunction and \
               not selectionFunction(childHwnd):
                continue
            results.append(childHwnd)
        return results

    return searchChildWindows(topHwnd)

def _sendNotifyMessage(hwnd, nofifyMessage):
    '''Send a notify message to a control.'''
    win32gui.SendMessage(win32gui.GetParent(hwnd),
                         win32con.WM_COMMAND,
                         _buildWinLong(nofifyMessage,
                                       win32api.GetWindowLong(hwnd,
                                                              win32con.GWL_ID)),
                         hwnd)

def _buildWinLong(high, low):
    '''Build a windows long parameter from high and low words.
    See http://support.microsoft.com/support/kb/articles/q189/1/70.asp
    '''
    # return ((high << 16) | low)
    return int(struct.unpack('>L', struct.pack('>2H', high, low)) [0])

def selectHoriz(hwnds, i):
    controls = []
    for x in hwnds:
        bbox = getBBox(x)
        controls.append([x]+bbox)
    x0 = sorted(list(set([x[1] for x in controls])))[i]
    return [x[0] for x in controls if x[1]==x0]

def selectVert(hwnds, i):
    controls = []
    for x in hwnds:
        bbox = getBBox(x)
        controls.append([x]+bbox)
    y0 = sorted(list(set([x[2] for x in controls])))[i]
    return [x[0] for x in controls if x[2]==y0]

def click(x,y):
    win32api.SetCursorPos((x,y))
    win32api.mouse_event(win32con.MOUSEEVENTF_LEFTDOWN,x,y,0,0)
    win32api.mouse_event(win32con.MOUSEEVENTF_LEFTUP,x,y,0,0)