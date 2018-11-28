import gzip
import sys

hap = sys.argv[1]
legend = sys.argv[2]
chrom = sys.argv[3]

hap_f = gzip.open(hap, 'rt', encoding='utf-8')
legend_f = gzip.open(legend, 'rt', encoding='utf-8')

legend_f.readline()  # header
snp  = legend_f.readline()
while snp != '':
    haps = hap_f.readline()
    sys.stdout.write('%s %s %s' % (chrom, snp.rstrip(), haps))
    snp = legend_f.readline()
