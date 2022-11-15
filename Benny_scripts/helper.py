# system imports
import os
# import sys

# third party imports

# local imports
from log_conf import log


generic_exception_string = "Please check input file and log file, something went wrong."


class StringNotFoundError(Exception):
    """Exception raised when a given string is not present in the target file"""
    def __init__(self, target_string, target_file):
        self.target_string = target_string
        self.target_file = target_file
        self.message = (
            f"It seems \"{target_string:s}\" was not present in the file\n\"{target_file:s}\"\n"
            "Check that the previous calculation didn't fail"
        )


def find_string_in_file(targetString, start_idx, end_idx, filedata, filePath):
    """wrapper that raises error if no substring can be found
    finds the first occurrence of a substring in memory mapped file
    """
    if targetString not in filedata:
        raise StringNotFoundError(targetString, filePath)

    location = filedata.find(targetString, start_idx+1, end_idx)

    if location != -1:
        log.debug(f"Found ({targetString=}) in {filePath}")
        return location
    else:
        raise StringNotFoundError(targetString, filePath)


def readfile(file_path, name=""):
    if not os.path.isfile(file_path):
        log.debug(f"The provided path is not valid:\n{file_path}\n")
        raise Exception(f"The provided path is not valid:\n{file_path}\n")

    try:
        with open(file_path, 'r') as fp:
            data = fp.read()
        return data

    except Exception as e:
        log.debug(f"Couldn't read the {name} file?\n{file_path}\n")
        raise e
