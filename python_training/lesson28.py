import re
namesRegex = re.compile(r'Agent \w+')
#namesRegex.findall('Agent Alice gave the secret documents to Agent Bob.')
#namesRegex.sub('REDACTED', 'Agent Alice gave the secret documents to Agent Bob.')

namesRegex = re.compile(r'Agent (\w)\w*')
namesRegex.findall('Agent Alice gave the secret documents to Agent Bob.')
namesRegex.sub(r'REDACTED \1****', 'Agent Alice gave the secret documents to Agent Bob.')
