""" helper.py"""

# system imports
import itertools as it
import os

# third party imports

# local imports


# -----------------------------------------------------------
# MEMORY MAPPED HELPER FUNCTIONS
# -----------------------------------------------------------

def readlines(mmFile, number_of_lines):
    """if we need to skip multiple lines"""
    for x in range(0, number_of_lines):
        mmFile.readline()


class StringNotFoundError(Exception):
    """Exception raised when a given string is not present in the target file"""
    def __init__(self, target_string, target_file):
        self.target_string = target_string
        self.target_file = target_file
        self.message = (
            f"It seems \"{target_string:s}\" was not present in the file\n\"{target_file:s}\"\n"
            "Check that the previous calculation didn't fail"
        )


def find_string_in_file(mmFile, filePath, targetString):
    """wrapper that raises error if no substring can be found
    finds the first occurrence of a substring in memory mapped file
    """
    location = mmFile.find(targetString.encode(encoding="utf-8"))

    if location == -1:
        raise StringNotFoundError(targetString, filePath)

    return location


def rfind_string_in_file(mmFile, filePath, targetString):
    """wrapper that raises error if no substring can be found
    finds the last occurrence of a substring in memory mapped file
    """
    location = mmFile.rfind(targetString.encode(encoding="utf-8"))

    if location == -1:
        raise StringNotFoundError(targetString, filePath)

    return location


def skip_back_lines(mm, numLines, startIndex):
    """gives the byte location numLines lines before
    the given byte location startIndex

    Factored out to simplify handling of n and offset"""

    for placeHolder in it.repeat(None, numLines):
        tempstartIndex = mm.rfind(b'\n', 0, startIndex)
        if tempstartIndex < 0:
            break
        startIndex = tempstartIndex

    return startIndex


def skip_forward_lines(mm, numLines, startIndex):
    """gives the byte location numLines lines after
    the given byte location startIndex

    Factored out to simplify handling of n and offset"""

    for placeHolder in it.repeat(None, numLines):
        tempstartIndex = mm.find(b'\n', startIndex + 1)
        if tempstartIndex == -1:
            break
        startIndex = tempstartIndex

    return startIndex

# -----------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------


def verify_file_exists(filePath):
    """x"""
    if not os.path.isfile(filePath):
        s = "The provided path {:s} does not appear to be a file, or it does not exist"
        raise FileNotFoundError(s.format(filePath))


# -----------------------------------------------------------
