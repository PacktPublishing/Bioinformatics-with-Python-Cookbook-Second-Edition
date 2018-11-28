import sys


sys.stdout.write('ID_1 ID_2 missing\n0 0 0 \n')
for line in sys.stdin:
    ind = line.rstrip()
    sys.stdout.write('%s %s 0\n' % (ind, ind))
    
