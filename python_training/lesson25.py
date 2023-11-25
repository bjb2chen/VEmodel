import re


batRegex = re.compile(r'Bat(wo)?man')
mo = batRegex.search('The Adventures of Batman')
mo.group()
'Batman'
mo = batRegex.search('The Adventures of Batwoman')
mo.group()
'Batwoman'
mo = batRegex.search('The Adventures of Batwowowowoman')
mo == None
True
phoneRegex = re.compile(r'\d\d\d-\d\d\d-\d\d\d\d')
mo = phoneRegex.search('My phone number is 415-555-1234. Call me tomorrow.')
mo.group()
'415-555-1234'
mo = phoneRegex.search('My phone number is 555-1234. Call me tomorrow.')
mo == None
True
phoneRegex = re.compile(r'(\d\d\d-)?\d\d\d-\d\d\d\d')
phoneRegex.search('My phone number is 415-555-1234. Call me tomorrow.')
<re.Match object; span=(19, 31), match='415-555-1234'>
phoneRegex.search('My phone number is 555-1234. Call me tomorrow.')
<re.Match object; span=(19, 27), match='555-1234'>


# asterisk * means zero or more times
batRegex = re.compile(r'Bat(wo)*man')
batRegex.search('The Adventures of Batman')
<re.Match object; span=(18, 24), match='Batman'>
batRegex.search('The Adventures of Batwoman')
<re.Match object; span=(18, 26), match='Batwoman'>
batRegex.search('The Adventures of Batowowowowowowowoman')
batRegex.search('The Adventures of Batwowowowowowowoman')
<re.Match object; span=(18, 38), match='Batwowowowowowowoman'>



batRegex = re.compile(r'Bat(wo)+man')
batRegex.search('THe Adventures of Batman')
batRegex.search('THe Adventures of Batman') == None
True
batRegex.search('THe Adventures of Batwoman') == None
False
batRegex.search('THe Adventures of Batwoman')
<re.Match object; span=(18, 26), match='Batwoman'>
batRegex.search('The Adventures of Batwowowowowowowoman')
<re.Match object; span=(18, 38), match='Batwowowowowowowoman'>


regex = re.compile(r'\+\*\?')
regex.search('I learned about +*? regex syntax')
<re.Match object; span=(16, 19), match='+*?'>
regex = re.compile(r'(\+\*\?)+')
regex.search('I learned about +*?+*?+*?+*?+*?+*? regex syntax')
<re.Match object; span=(16, 34), match='+*?+*?+*?+*?+*?+*?'>


haRegex = re.compile(r'(Ha){3}')
haRegex.search('He said "HaHaHa"')
<re.Match object; span=(9, 15), match='HaHaHa'>
phoneRegex = re.compile(r'((\d\d\d-)?\d\d\d-\d\d\d\d(,)?){3}')
phoneRegex.search('My numbers are 415-555-1234, 555-4242, 212-555-0000')
phoneRegex.search('My numbers are 415-555-1234,555-4242,212-555-0000')
<re.Match object; span=(15, 49), match='415-555-1234,555-4242,212-555-0000'>



haRegex = re.compile(r'(Ha){3,5}')
haRegex.search('He said "HaHaHa"')
<re.Match object; span=(9, 15), match='HaHaHa'>
haRegex.search('He said "HaHaHaHaHa"')
<re.Match object; span=(9, 19), match='HaHaHaHaHa'>
haRegex = re.compile(r'(Ha){3,}'
                     )
digitRegex = re.compile(r'(\d){3,5}')
digitRegex.search('1234567890')
<re.Match object; span=(0, 5), match='12345'>
# notice how python matched 5 instead of 3, is because python is greedy matcher
digitRegex = re.compile(r'(\d){3,5}?')
digitRegex.search('1234567890')
<re.Match object; span=(0, 3), match='123'>
# the ? enforces non-greedy lazy match

# recap: ? says zero or one times, * zero or more times, + one or more times
