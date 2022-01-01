# This file is generated by wxPython's SIP generator.  Do not edit by hand.
#
# Copyright: (c) 2020 by Total Control Software
# License:   wxWindows License

"""
The ``wx.media`` module provides a widget class that allows displaying various
types of media, such as video and audio files and streaming, using native
system components.  The wxWidgets media classes are an optional part of the
build so it may not always be available on your build of wxPython.
"""

from ._media import *

import wx

EVT_MEDIA_LOADED = wx.PyEventBinder( wxEVT_MEDIA_LOADED )
EVT_MEDIA_STOP = wx.PyEventBinder( wxEVT_MEDIA_STOP )
EVT_MEDIA_FINISHED = wx.PyEventBinder( wxEVT_MEDIA_FINISHED )
EVT_MEDIA_STATECHANGED = wx.PyEventBinder( wxEVT_MEDIA_STATECHANGED )
EVT_MEDIA_PLAY = wx.PyEventBinder( wxEVT_MEDIA_PLAY )
EVT_MEDIA_PAUSE = wx.PyEventBinder( wxEVT_MEDIA_PAUSE )

MEDIABACKEND_DIRECTSHOW = "wxAMMediaBackend"
MEDIABACKEND_MCI        = "wxMCIMediaBackend"
MEDIABACKEND_QUICKTIME  = "wxQTMediaBackend"
MEDIABACKEND_GSTREAMER  = "wxGStreamerMediaBackend"
MEDIABACKEND_REALPLAYER = "wxRealPlayerMediaBackend"
MEDIABACKEND_WMP10      = "wxWMP10MediaBackend"

