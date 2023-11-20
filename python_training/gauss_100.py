total = 0
for num in range(101):
    total = total + num
print(total)

# the way Gauss figured it out was to add 1 + 99, 2 + 98, ...
# for fifty pairs, so 50 * 100 (the value of the sum of each pair)
# + a final 50 for the middle 50 which is not with a pair
# you get = 5050.
