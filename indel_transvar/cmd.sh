#/lustre/work/user/chenmj/software/Anaconda3-2021.05/bin/transvar ganno -l input.txt --reference /opt/anaconda2/lib/python2.7/site-packages/transvar/transvar.download/hg19.fa --ensembl /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ensembl.gtf.gz.transvardb --ucsc /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ucsc.txt.gz.transvardb --gseq > input.txt.transvar
rsync -aivHLP chen.minjun@192.168.11.254:/lustre/work/user/chenmj/analysis/transvar/input.txt.transvar .