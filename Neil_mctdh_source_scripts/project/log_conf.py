import logging

# -----------------------------------------------------------
# LOGGING PREPERATIONS
# -----------------------------------------------------------
# predefined levels for logging
# CRITICAL 50
# ERROR    40
# WARNING  30
# logging.FLOW = 25
# INFO     20
# DEBUG    10
# NOTSET   0
# -----------------------------------------------------------

# add names
# logging.addLevelName(logging.FLOW, "FLOW")


class MyLogger(logging.Logger):
    def flow(self, message, *args, **kwargs):
        if self.isEnabledFor(logging.FLOW):
            self._log(logging.FLOW, message, args, **kwargs)


logging.setLoggerClass(MyLogger)
log = logging.getLogger(__name__)

logging.basicConfig(
    format="%(asctime)-13s [%(levelname)s] %(funcName)s: %(message)s",
    # datefmt='%m/%d/%Y %I:%M:%S %p',
    datefmt='%d %I:%M:%S ',
    # level=logging.FLOW,
    level=logging.INFO,
    # level=logging.DEBUG,
    # level=logging.LOCK,
)


def setLevelCritical():
    log.setLevel(logging.CRITICAL)

def setLevelError():
    log.setLevel(logging.ERROR)

def setLevelWarning():
    log.setLevel(logging.WARNING)

def setLevelFlow():
    log.setLevel(logging.FLOW)

def setLevelInfo():
    log.setLevel(logging.INFO)

def setLevelDebug():
    log.setLevel(logging.DEBUG)