"""
timer.py
Timer class for timing nested sections of code
"""

# python libraries
import os
import sys
import traceback
import time

# GLOBALS
_RASMUS_TIMER = None
_GLOBAL_NOTES = None


#=============================================================================

class Timer(object):
    """Timer class"""

    def __init__(self, stream=sys.stderr, maxdepth=1e1000):
        self.reset()
        self.streams = [(stream, maxdepth)]
        self.show_errors = True
        self.show_warnings = True
        self.quiets = 0


    def start(self, msg=""):
        """Start a new timer"""

        if msg != "":
            self.indent()
            self._write("BEGIN %s:\n" % msg)
        self.msg.append(msg)
        self.flush()
        self.starts.append(time.time())


    def time(self):
        """Get current duration of the timer"""

        return self.starts[-1] - time.clock()


    def stop(self):
        """Stop last created timer and return duration in seconds"""

        duration = time.time() - self.starts.pop()
        msg = self.msg.pop()
        if msg != "":
            self.indent()

            if duration > 3600:
                pretty = "%.1fh" % (duration / 3600.)
            elif duration > 60:
                pretty = "%.1fm" % (duration / 60.)
            else:
                pretty = "%.3fs" % duration

            if duration > .1:
                secs = "%.3fs" % duration
            else:
                secs = "%.3es" % duration

            self.write("END   %s: %s (%s)\n" % (msg, pretty, secs))
        self.flush()
        return duration


    def log(self, *text):
        """Write message to timer stream
        Message will be written with current indentation level"""

        self.indent()
        for i in text:
            self._write("%s " % str(i))
        self._write("\n")
        self.flush()


    def log_exact(self, text):
        """Write extact string 'text' to timer output stream
        with no additional indentation"""

        self._write(text)
        self.flush()


    def warn(self, text, offset=0):
        """Write warning message to timer output stream"""

        filename, lineno, _, _ = traceback.extract_stack()[-2-offset][:2]
        filename = os.path.basename(filename)

        if self.show_warnings:
            self.indent()
            self._write("WARNING: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()


    def error(self, text, offset=0):
        """Write error message to timer output stream"""

        filename, lineno, _, _ = traceback.extract_stack()[-2-offset]
        filename = os.path.basename(filename)

        if self.show_errors:
            self.indent()
            self._write("ERROR: %s, line %d: %s\n" % (filename, lineno, text))
            self.flush()


    def indent(self):
        """Write current indentation level to timer output stream"""
        for _ in range(self.depth()):
            self._write("  ")


    def reset(self):
        """Stop all timers"""
        self.msg = []
        self.starts = []


    def depth(self):
        """Get current number of running timers"""
        return len(self.msg)


    def _write(self, text):
        """Private function for writing to output stream"""
        for stream, maxdepth in self.streams:
            if self.depth() < maxdepth and \
               self.quiets == 0:
                stream.write(text)


    def write(self, text):
        """Write text to output stream"""
        self._write(text)
        self.flush()


    def flush(self):
        """Flush all timers"""
        for stream, _ in self.streams: # stream, maxdepth
            stream.flush()


    def add_stream(self, stream, maxdepth=1e1000):
        """Add stream with maxdepth"""
        self.streams.append((stream, maxdepth))


    def remove_stream(self, stream):
        """Remove stream"""
        self.streams = [x for x in self.streams if x[0] != stream]


    def suppress(self):
        """Calling this function will suppress timer output messages until
           unsuppress() is called.

           If suppress() is called multiple times,  unsuppress() must be called
           an equal number of times to resume timer  output.  This is useful for
           nesting suppress/unsuppress."""
        self.quiets += 1


    def unsuppress(self):
        """Calling this function will resume timer output messages that were
           disabled with suppress().

           If suppress() is called multiple times,  unsuppress() must be called
           an equal number of times to resume timer  output.  This is useful for
           nesting suppress/unsuppress."""
        self.quiets = max(self.quiets - 1, 0)


#=============================================================================
# Convenience functions

def global_timer():
    """Global timer"""
    global _RASMUS_TIMER
    if _RASMUS_TIMER is None:
        _RASMUS_TIMER = Timer()
    return _RASMUS_TIMER


def log(*text):
    """Write text to global timer"""
    return global_timer().log(*text)
logger = log


def log_exact(text):
    """Write exact text to global timer"""
    return global_timer().log_exact(text)


def tic(msg=""):
    """Start global timer"""
    return global_timer().start(msg)


def toc():
    """Stop global timer"""
    return global_timer().stop()


def indent():
    """Indent global timer"""
    return global_timer().indent()


def warn(text, offset=0):
    """Write warning message to global timer"""
    return global_timer().warn(text, offset+1)


def error(text, offset=0):
    """Write error message to global timer"""
    return global_timer().error(text, offset+1)


def note(*text):
    """Write text to global note file"""
    print >>notefile(), " ".join(text)


def noteflush():
    """Flush global note file"""
    return notefile().flush()


def notefile(out=None):
    """Global note file"""
    global _GLOBAL_NOTES

    if out is None:
        out = file("/dev/null", "w")
    if _GLOBAL_NOTES is None:
        _GLOBAL_NOTES = out
    return _GLOBAL_NOTES


#=============================================================================
# debugging info

def current_file(offset=0, abbrv=True):
    """Get traceback filename"""
    filename = traceback.extract_stack()[-2-offset][0]
    if abbrv:
        filename = os.path.basename(filename)
    return filename


def current_line(offset=0):
    """Get traceback line number"""
    lineno = traceback.extract_stack()[-2-offset][1]
    return lineno


def current_func(offset=0):
    """Get traceback function"""
    func = traceback.extract_stack()[-2-offset][2]
    return func


def current_code(offset=0):
    """Get traceback code"""
    code = traceback.extract_stack()[-2-offset][3]
    return code
