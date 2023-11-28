Python 3.11.5 (tags/v3.11.5:cce6ba9, Aug 24 2023, 14:38:34) [MSC v.1936 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
import re
beginsWithHelloRegex = re.compile(r'^Hello')
beginsWithHelloRegex.search('Hello there!')
<re.Match object; span=(0, 5), match='Hello'>
beginsWithHelloRegex.search('He said "Hello!"')
beginsWithHelloRegex.search('He said "Hello!"') == None
True
endsWithWorldRegex = re.compile(r'world!$')
endsWithWorldRegex.search('Hello world!')
<re.Match object; span=(6, 12), match='world!'>
endsWithWorldRegex.search('Hello world!srghsdljkahfdjash')
endsWithWorldRegex.search('Hello world!srghsdljkahfdjash') == None
True


allDigitsRegex = re.compile(r'^\d+$')
allDigitsRegex.search('2343729479832754')
<re.Match object; span=(0, 16), match='2343729479832754'>
allDigitsRegex.search('2343729479x32754') == None
True
# ^both$ means that the entire string must match
atRegex = re.compile(r'.at')
atRegex.findall('The cat in the hat sat on the flat mat.')
['cat', 'hat', 'sat', 'lat', 'mat']


atRegex = re.compile(r'.{1,2}at')
atRegex.findall('The cat in the hat sat on the flat mat.')
[' cat', ' hat', ' sat', 'flat', ' mat']
'First Name: Al Last Name: Sweigart'
'First Name: Al Last Name: Sweigart'
nameRegex = re.compile(r'First Name: (.*) Last Name: (.*)')
nameRegex.findall('First Name: Al Last Name: Sweigart')
[('Al', 'Sweigart')]
serve = '<To serve humans> for dinner.>'
nongreedy = re.compile(r'<(.*?)>')
nongreedy.findall(serve)
['To serve humans']
greedy = re.compile(r'<(.*)>')
greedy.findall(serve)
['To serve humans> for dinner.']
prime = 'Serve the public trrust. \nProtect the innocent.\nUpload the law.'
print(prime)
Serve the public trrust. 
Protect the innocent.
Upload the law.
dotStar = re.compile(r'.*')
dotStar.search(prime)
<re.Match object; span=(0, 25), match='Serve the public trrust. '>
dotStar = re.compile(r'.*', re.DOTALL)
dotStar.search(prime)
<re.Match object; span=(0, 63), match='Serve the public trrust. \nProtect the innocent.\>
vowelRegex = re.compile(r'[aeiou]')
vowelRegex.search('Al, why does your programming booko talk about RoboCop so much?')
<re.Match object; span=(9, 10), match='o'>
vowelRegex.findall('Al, why does your programming booko talk about RoboCop so much?')
['o', 'e', 'o', 'u', 'o', 'a', 'i', 'o', 'o', 'o', 'a', 'a', 'o', 'u', 'o', 'o', 'o', 'o', 'u']
vowelRegex = re.compile(r'[aeiou]', re.I)
vowelRegex.findall('Al, why does your programming booko talk about RoboCop so much?')
['A', 'o', 'e', 'o', 'u', 'o', 'a', 'i', 'o', 'o', 'o', 'a', 'a', 'o', 'u', 'o', 'o', 'o', 'o', 'u']
