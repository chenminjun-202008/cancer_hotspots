# 下载数据库
`wget https://www.cancerhotspots.org/files/hotspots_v2.xls`
`wget http://download.cbioportal.org/cancerhotspots/cancerhotspots.v2.maf.gz`

---
# 合并snv/indel到1个文件，方便从maf一起提取
`Rscript scripts/hotspots_combine.R hotspots_v2.xls hotspots_v2.txt`

# 将maf中绝对坐标补充到hotspots_v2
`/mnt/titan01/Orca/cupcake/softwares/python-3.9/bin/python3 scripts/get_cancerhotspots.py -m cancerhotspots.v2.maf -s hotspots_v2.txt -o hotspots_v2.add.txt`

# maf有部分位点不存在，使用transvar进行补充


```
rsync -aivHLP hotspots_v2.add.txt.err_for_transvar chen.minjun@192.168.11.254:/lustre/work/user/chenmj/analysis/transvar
/lustre/work/user/chenmj/software/Anaconda3-2021.05/bin/transvar panno -l hotspots_v2.add.txt.err_for_transvar -t 1 -m 2 --reference /opt/anaconda2/lib/python2.7/site-packages/transvar/transvar.download/hg19.fa --ensembl /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ensembl.gtf.gz.transvardb > transvar-ucsc.err.txt
less transvar-ucsc.err.txt |awk '$2=="." {print $1}'|sed 's/|/\t/g' > tmp2
/lustre/work/user/chenmj/software/Anaconda3-2021.05/bin/transvar panno -l tmp2 -m 2 --reference /opt/anaconda2/lib/python2.7/site-packages/transvar/transvar.download/hg19.fa --ensembl /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ensembl.gtf.gz.transvardb --ucsc /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ucsc.txt.gz.transvardb >> transvar-ucsc.err.txt

rsync -aivHLP chen.minjun@192.168.11.254:/lustre/work/user/chenmj/analysis/transvar/transvar-ucsc.err.txt ./
```

#补充位点并合并到一起
```
Rscript scripts/add_err_site.R  hotspots_v2.add.txt.err transvar-ucsc.err.txt transvar-ucsc.err.transvar.txt

awk 'NR==1||FNR>1' hotspots_v2.add.txt transvar-ucsc.err.transvar.txt > hotspots_v2.addAll.txt
```

---
# indel位点，使用transvar转化为非“-”形式

```bash
less hotspots_v2.addAll.txt |awk '$4=="-" {OFS="\t"; print "chr"$1":"$2"_"$3"ins"$5,$1,$2,$3,$4,$5}' > indel_transvar/input.txt
less hotspots_v2.addAll.txt |awk '$5=="-" {OFS="\t";print "chr"$1":"$2"_"$3"del"$4,$1,$2,$3,$4,$5}' >> indel_transvar/input.txt


/lustre/work/user/chenmj/software/Anaconda3-2021.05/bin/transvar ganno -l input.txt --reference /opt/anaconda2/lib/python2.7/site-packages/transvar/transvar.download/hg19.fa --ensembl /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ensembl.gtf.gz.transvardb --ucsc /lustre/work/user/chenmj/software/Anaconda3-2021.05/lib/python3.8/site-packages/transvar/transvar.download/hg19.ucsc.txt.gz.transvardb --gseq > input.txt.transvar
rsync -aivHLP chen.minjun@192.168.11.254:/lustre/work/user/chenmj/analysis/transvar/input.txt.transvar .


Rscript scripts/get_indel_refAlt.R hotspots_v2.addAll.txt indel_transvar/input.txt indel_transvar/input.txt.transvar hotspots_v2.final.txt
```
