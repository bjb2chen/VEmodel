import re
namesRegex = re.compile(r'Agent \w+')
#namesRegex.findall('Agent Alice gave the secret documents to Agent Bob.')
#namesRegex.sub('REDACTED', 'Agent Alice gave the secret documents to Agent Bob.')

namesRegex = re.compile(r'Agent (\w)\w*')
namesRegex.findall('Agent Alice gave the secret documents to Agent Bob.')
namesRegex.sub(r'REDACTED \1*\*\*\*', 'Agent Alice gave the secret documents to Agent Bob.')


re.compile(r'''
	(\d\d\d-)|   # area code
	(\(\d\d\d\) )# -or- area code with parens and no dash
	-        # first dash
	\d\d\d   # first 3 digits
	-        # second dash
	\d\d\d\d
	\sx\d{2,4} # etension, like x1234
	''', re.VERBOSE | re.DOTALL | re.IGNORECASE)
