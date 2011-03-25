#!/usr/bin/env python
#
# Coalesces a list such as "1 2 3 4 5 10 15 19 20 21" to "1-5; 10; 15; 19-21"
# for use with the print_res keyword of MMPBSA.py.
import sys

# Read and parse numbers
data = " ".join(sys.stdin.readlines()).split()
numbers = [int(x) for x in data]

def print_range(s, e):
    if s != e:
        print "%d-%d;" % (s, e),
    else:
        print "%d;" % s,
        
start = numbers[0]
end = numbers[0]
last = start
for i in xrange(1, len(numbers)):
    num = numbers[i]
    if num == last + 1:
        end = num
    else:
        print_range(start, end)
        start = end = num
    last = num
print_range(start, end)