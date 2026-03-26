zcat  /EClab/Project/Snow_leopard/Future/06.Know_new/00.VCF/All.70.screen.vcf.gz  |grep -m 1 CHRO |tr '\t' '\n' |tail -n +10 > VCF.order.list
perl /EClab/Project/Snow_leopard/Future/06.Know_new/INFO/Final/step2.ad.info.pl  /EClab/Project/Snow_leopard/Future/06.Know_new/INFO/Final/All.176.info.new.txt  VCF.order.list VCF.order.ad.txt
awk -F '\t' '{print $1"\t"$2}' VCF.order.ad.txt > Dsuite.cluster
