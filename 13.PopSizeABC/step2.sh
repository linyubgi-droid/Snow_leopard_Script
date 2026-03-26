mkdir -p /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS
for i in {1..18} ;do echo zcat /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS.screen.v41.change.vcf.gz \| awk \'\$1\~\/\^\#\/ \|\| \$1\~\/chr$i\/ \{print \$0\}\'  \| gzip \-f \> chr$i".QS.vcf.gz" ;done> /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/QS/step1.split.sh
mkdir -p /EClab/Project/Snow_leopard/Future/06.Know_new/16.PopsizeABC_69/QS/res_final
