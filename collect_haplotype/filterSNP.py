import sys
from numpy import random as rd


def readfile(sfile):

	lfilter = []

	ntotal = 0
	with open(sfile, 'r') as f:
		for l in f:
			l = l.strip()
			if l.startswith('#'):
				lfilter.append(l + '\n')
				continue
			
			r = l.split('\t')
			if len(r[3]) != 1 or len(r[4]) != 1:
				continue
			lfilter.append(l + '\n')
								
			ntotal += 1

	print 'total # of snps in chr10:', ntotal


        f = open('chr10_1kgenomes_SNPs.vcf', 'w')
        for l in lfilter:
                f.write(l)

        f.close()




			


if __name__ == '__main__':
	readfile('/data/diybu/1000genomes/chr_vcf/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf')
	


