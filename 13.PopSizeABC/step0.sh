cp /EClab/Project/Snow_leopard/Future/06.Know_new/00.VCF/01.region_vcf_69/bin/QS.list list_indiv_QS.txt
ln -sf  /EClab/Project/Snow_leopard/Future/06.Know_new/00.VCF/01.region_vcf_69/QS/QS.screen.vcf.gz /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS
perl /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/bin/vcf4.2to4.1.pl /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS.screen.vcf.gz /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS.screen.v41.vcf.gz
zcat /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS.screen.v41.vcf.gz | sed 's/0|0/0\/0/g; s/0|1/0\/1/g; s/1|1/1\/1/g'  |gzip -f > /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS.screen.v41.change.vcf.gz
cat list_indiv_QS.txt | awk '{print "QS",$0,"0 0 1 -999"}' >  QS.indiv.ped
