
import sys
import os
import argparse
import random
import time
import pickle as pk
import copy
import numpy as np
import scipy.stats as st#import chisquare
import matplotlib.pyplot as plt
#import cPickle as pickle


def readfile(sfile):

        dsnpFreq = {} #whole DB snp freq
	dsnpFreqSub = {} #subset DB snp freq
	dsnpForm = {} #snp ref alt allele
        #print 'lmark', lmark
	fout = open('chr10Subset.vcf', 'w')
        with open(sfile, 'r') as f:
                for l in f:
                        l = l.strip()
                        if l.startswith('##'):
                                continue

                        if l.startswith('#CHROM'):
				r = l.split()
				N = len(r)
				random.seed(716)
				lsubIdx = sorted(random.sample(range(9, N), 500)) #500 is the size of subset DB	
				lsubid = '\t'.join(r[:9] + [r[i] for i in lsubIdx])
				fout.write(lsubid + '\n')	
                                continue


                        r = l.split()
			srsid = ' '.join(r[:3])
                        #srsid = r[2]
			#if srsid == '.':
			#	continue
			#elif 'rs' not in srsid and 'esv' not in srsid:
			#	print l
			if len(r[3]) != 1 or len(r[4]) != 1:
                                continue 

			dsnpForm[srsid] = r[3:5]
			

                       # linfo = r[7].split(';')


                       # for i, elem in enumerate(linfo):
                       #         if 'AF=' in elem:
                       #                 faf = float(elem.split('=')[1])
			
			cnt = 0
			cnt_sub = 0
			for idx in range(9, len(r)):
				if '1' in r[idx]:
					#cnt += 1
					ntmp = r[idx].count('1')
					cnt += ntmp
					if idx in lsubIdx:
						#cnt_sub += 1
						cnt_sub += ntmp
			faf = float(cnt)/(2*len(r[9:]))
			faf_sub = float(cnt_sub)/(2*len(lsubIdx))


			#dsnpFreq[srsid] = cnt
			dsnpFreq[srsid] = faf

			#faf_sub = ' '.join(r[9:]).count('1')/(2*float(len(lsubIdx)))
			#dsnpFreqSub[srsid] = cnt_sub
			dsnpFreqSub[srsid] = faf_sub
		
			lnew = '\t'.join(r[:7] + ['AF='+str(faf_sub)+';VT=SNP'] + [r[8]] + [r[i] for i in lsubIdx])
			fout.write(lnew + '\n')



	fout.close()

	dsnpFreq['TOTAL'] = 2*len(r[9:])
	dsnpFreqSub['TOTAL'] = 2*len(lsubIdx)

        output = open('dsnpFreqWholeDB.pk', 'wb')
	pk.dump(dsnpFreq, output)
	output.close()

	output_sub = open('dsnpFreqSubset.pk', 'wb')
	pk.dump(dsnpFreqSub, output_sub)
	output_sub.close()


	output_allele = open('dsnpAllele.pk', 'wb')
	pk.dump(dsnpForm, output_allele)
	output_allele.close()


if __name__ == '__main__':


	path = '/data2/diybu/beaconHaplopAttack/haploblock/'
	dtype = readfile('/data/diybu/1000genomes/chr_vcf/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf')

