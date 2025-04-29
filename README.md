# ATAC


```

```

```

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
zcat refFlat.txt.gz | perl -alne '{next if /^#/;if($F[3] eq "+"){$start=$F[4]-2500;$end=$F[4]+2500}else{$start=$F[5]-2500;$end=$F[5]+2500}print join("\t",$F[2],$start,$end,$F[12],0,$F[3])}' | sort -u > hg38.refseq.tss.bed


 wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz^C
 zcat refFlat.txt.gz | perl -alne '{next if /^#/;if($F[3] eq "+"){$start=$F[4]-2500;$end=$F[4]+2500}else{$start=$F[5]-2500;$end=$F[5]+2500}print join("\t",$F[2],$start,$end,$F[12],0,$F[3])}' | sort -u > hg19.refseq.tss.bed
```

```
python run_atac_visualization.py -b . -o Visualization_Output  --tss_bed hg19.tss.bed  --homer_genome hg19
```

```
python run_atac_qc.py -b . -o QC_Results --tss_bed hg19.tss.bed --cores 8 --parafly_jobs 4
```
