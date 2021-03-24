

for ((i=20000;i<=135700000-150000;i=i+150000));
do
	#echo "10:$(($i - 20000))-$(($i + 150000))"
	perl vcf_to_ped_converter.pl -vcf chr10_1kgenomes_SNPs.vcf -sample_panel_file /data/diybu/1000genomes/integrated_call_samples_v3.20130502.ALL.panel -region "10:$(($i - 20000))-$(($i + 150000))" -population GBR -population FIN -population ACB -population ASW -population BEB -population CDX -population CEU -population CHB -population CHS -population CLM -population ESN -population GIH -population GWD -population IBS -population ITU -population JPT -population KHV -population LWK -population MSL -population MXL -population PEL -population PJL -population PUR -population STU -population TSI -population YRI 


done
