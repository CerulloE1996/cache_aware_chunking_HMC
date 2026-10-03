## Prints the X session's input idle time in seconds (DISPLAY :1), via libXss.
import ctypes, ctypes.util
x11 = ctypes.cdll.LoadLibrary(ctypes.util.find_library("X11"))
xss = ctypes.cdll.LoadLibrary(ctypes.util.find_library("Xss"))
class XScreenSaverInfo(ctypes.Structure):
    _fields_ = [("window", ctypes.c_ulong), ("state", ctypes.c_int), ("kind", ctypes.c_int), ("til_or_since", ctypes.c_ulong), ("idle", ctypes.c_ulong), ("eventMask", ctypes.c_ulong)]
x11.XOpenDisplay.restype = ctypes.c_void_p
x11.XDefaultRootWindow.restype = ctypes.c_ulong
x11.XDefaultRootWindow.argtypes = [ctypes.c_void_p]
xss.XScreenSaverAllocInfo.restype = ctypes.POINTER(XScreenSaverInfo)
xss.XScreenSaverQueryInfo.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.POINTER(XScreenSaverInfo)]
display = x11.XOpenDisplay(b":1")
info = xss.XScreenSaverAllocInfo()
xss.XScreenSaverQueryInfo(display, x11.XDefaultRootWindow(display), info)
print(int(info.contents.idle / 1000))
