spam = 'Hello world!'

print(spam.upper())
print(spam.lower())

# answer = input()
# print('The answer you inputted is: ' + answer)

# if answer == 'yes':
#     print('Playing again')

str1 = 'Hello'
print(str1.isupper())
print(str1.istitle())

','.join(['cats', 'rats', 'bats'])
print(','.join(['cats', 'rats', 'bats']))

print('Hello'.rjust(10))
print('Hello'.ljust(20))
print('Hello'.rjust(20, '*'))
print('Hello'.center(20,'='))
print('Hello'.center(20))
name = 'Al'
print(name.center(20, '='))

'SpamSpamBaconSpamEggSpamSpam'.strip('ampS')
print('SpamSpamBaconSpamEggSpamSpam'.strip('ampS')) # this will strip all a, m, p, S chars in the str until it reaches something not so

spam = 'Hello there!'
spam = spam.replace('e', 'XYZ')
print(spam)

import pyperclip
