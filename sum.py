#!/usr/bin/env python

from sys import stdin, stderr

lines = stdin.readlines()

total = 0
num_values = 0
for line in lines:
    try:
        n = float(line.strip())
        total += n
        num_values += 1
    except ValueError:
        continue

print(f"Summed {num_values} values in {len(lines)} lines of input.", file=stderr)
print(total)
