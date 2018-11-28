from __future__ import print_function
import sys

stats = sys.argv[1] == '--stats=yes'
input_file = sys.argv[2]
output_file = sys.argv[3]

print(sys.argv)

num_recs = 0
positions = set()
with open(input_file) as f:
    f.readline()
    for l in f:
        num_recs += 1
        toks = l.rstrip().split('\t')
        positions.add(int(toks[1]))
        positions.add(int(toks[2]))


with open(output_file, 'w') as w:
    w.write('Records: %d\n' % num_recs)
    if stats:
        w.write('Minimum position: %d\n' % min(positions))
        w.write('Maximum position: %d\n' % max(positions))
