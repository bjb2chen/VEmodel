# system imports
import os


# third party imports


# local imports


# -----------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------

def verify_file_exists(filePath):
    """x"""
    if not os.path.isfile(filePath):
        s = "The provided path {:s} does not appear to be a file, or it does not exist"
        raise FileNotFoundError(s.format(filePath))
