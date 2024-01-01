#! python3

import re
import shutil


# TODO: Create a regex object for phone numbers
re.compile(r'''

# 415-555-0000, 555-0000, (415) 555-0000, 555-0000 ext. 12345, x12345

(\d\d\d\)|(\(\d\d\d\)))?		# area code (optiional)
		# first separator
		# first 3 digits
		# separator
		# last 4 digits
		# extension (optional)

''', re.VERBOSE)

# TODO: Create a regex object for email addresses

# TODO: Get the text off the clipboard

# TODO: Extract the email/phone from this text

# TODO: Copy the extracted phone/email to the clipboard