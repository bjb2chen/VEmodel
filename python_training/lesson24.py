import re
phoneNumRegex = re.compile(r'\d\d\d-\d\d\d-\d\d\d\d')
phoneNumRegex.search('My number is 415-555-4242')
<re.Match object; span=(13, 25), match='415-555-4242'>
mo = phoneNumRegex.search('My number is 415-555-4242')
mo.group()
'415-555-4242'
phoneNumRegex = re.compile(r'(\d\d\d)-(\d\d\d-\d\d\d\d)')
mo = phoneNumRegex.search('My number is 415-555-4242')
mo.group()
'415-555-4242'
mo.group(1)
'415'
mo.group(2)
'555-4242'
phoneNumRegex = re.compile(r'\(\d\d\d\) \d\d\d-\d\d\d\d')
mo = phoneNumRegex.search('My number is (415) 555-4242')
mo.group()
'(415) 555-4242'


batRegex = re.compile(r'Bat(man|mobile|copter|bat)')
mo = batRegex.search('Batmobile lost a wheel')
mo.group()
'Batmobile'
mo = batRegex.search('Batmotorcycle lost a wheel')
mo == None
True
mo = batRegex.search('Batmobile lost a wheel')
mo.group(1)
'mobile'
mo.group()
'Batmobile'
