import logging
from PyQt5.QtGui import QColor
from PyQt5.QtCore import QObject, pyqtSignal

class QTextEditLogger(logging.Handler, QObject):
    append_signal = pyqtSignal(str, str)

    def __init__(self, widget):
        logging.Handler.__init__(self)
        QObject.__init__(self)
        self.widget = widget
        self.widget.setReadOnly(True)
        self.append_signal.connect(self.append_text)

    def emit(self, record):
        msg = self.format(record)
        if record.levelno == logging.DEBUG:
            color = 'gray'
        elif record.levelno == logging.INFO:
            color = 'blue'
        elif record.levelno == logging.WARNING:
            color = 'yellow'
        elif record.levelno == logging.ERROR:
            color = 'red'
        elif record.levelno == logging.CRITICAL:
            color = 'darkred'
        # Emit signal with message and color
        self.append_signal.emit(msg, color)

    def append_text(self, msg, color):
            self.widget.setTextColor(QColor(color))
            self.widget.append(msg)
            self.widget.setTextColor(QColor('black'))