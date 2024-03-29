HTSLIB=/usr/local/genome/htslib-1.19.1

all:
	g++ -g2 -Wall -I$(HTSLIB) -o mtdoc mtmain.cc  -L$(HTSLIB) -lhts -lpthread
mtest:
	export LD_LIBRARY_PATH=$(HTSLIB); ./mtdoc --bam /ssd3/parabricks/Sample_NIST7035.bam --bed Twist_Comprehensive_Exome_Covered_Targets_hg38.xhmm.interval_list --out output

