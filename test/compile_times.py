import sys
from parse import parse

with open(sys.argv[1]) as f:
    lines = f.readlines()
desired_lines = lines[2:len(lines):4]

values = []

for line in desired_lines:
    values.append(int(parse('Algorithm took {} milliseconds\n', line)[0]))

print(values)