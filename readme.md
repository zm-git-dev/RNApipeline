###This is a RNA Pipeline box, include three different RNA pipeline

* `Author:` yueyao@genomics.cn
* `RNAseq`
* `RNAref`
* `RNAdenovo`

RNA流程搭建主要包括RNAseq，RNAref，RNAdenovo三个部分，流程主要功能已经实现，可以针对动物和植物进行不同的分析

* RNAseq包括过滤，比对到基因组，表达量的计算，差异表达分析，GO富集分析，KEGG富集分析，共表达网络计算(样品数大于6时)
* RNAref包括过滤，比对到基因组，新转录本的预测，可变剪切分析，SNP和INDEL分析，表达量的计算，差异表达分析，转录因子预测，蛋白互作网络分析，基因融合分析，CIRCOS图，植物抗病基因注释或者真菌致病基因注释，聚类分析，GO富集分析，KEGG富集分析
* RNAdenovo包括转录本的拼接，SSR分析，CDS预测，基因的定量，SNP分析，基因的功能注释，差异表达分析，GO富集分析，KEGG富集分析，转录因子预测，植物抗病基因注释，蛋白互作网络分析
* 这里面包括集群的版本和本地的版本，本地如果配置好软件和数据库以及运行所需要的脚步是可以直接运行的
