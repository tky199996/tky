#qiime2分析流程
#文件准备:
manifest.txt(导入测序数据)   samplemetadata.txt(后续分组比较所需文件)
#下载qiime2:两个版本都需要:
qiime2-2023.9-amplicon与qiime2-2023.9-shotgun
#数据导入
time qiime tools import \
--type 'SampleData[SequnecesWithQuality]' \
--input-path manifest.txt \
--output-path single-end-demux.qza \
--input-format SingleEndFastqManifestPhred33V2 
##DADA2降噪
#创建文件夹
mkdir dada 
qiime dada2 denoise-single \
--i-demultiplexed-seqs Single-end-demux.qza \
--p-trim-left 0 \
--p-trunc-len 0 \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza \
--o-denoising-stats stats.qza 

#特征表可视化
qiime featue-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file samplemetadata.txt 

#创建ASV表，构建系统进化树
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \ （代表序列）
--o-alignment aligned-rep-seqs.qza \  （序列对齐）
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#α多样性与β多样性分析
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 40000 #根据table.qzv文件确定抽平丰度，一般可以选择最小丰度值的90%-95% \
--m-metadata-file samplemetadata.txt \
--output-dir core-metrics-results

#Alpha 多样性指数的组间差异分析
# faith-pd
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
--m-metadata-file samplemetadata.txt \
--o-visualization core-metrics-results/faith-pd-group-significance.qzv
相同组织，不同基因型间faith-pd指数无差异；不同组织间faith-pd指数差异显著
# evenness
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/evenness_vector.qza \
--m-metadata-file samplemetadata.txt \
--o-visualization core-metrics-results/evenness-group-significance.qzv
相同组织，不同基因型间faith-pd指数无差异；不同组织间faith-pd指数差异显著

#Beta多样性指数的组间差异分析
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file sample-metadata.txt \
--m-metadata-column group \
--o-visualization core-metrics-results/unweighted-unifrac-genetype-position-significance.qzv \
--p-pairwise
十二指肠不同基因型间qvalue < 0.05，盲肠不同基因型间beta多样性指数无差异。

#Alpha稀释曲线
qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 56,768 \#最大采样深度，一般根据table.qzv取frequency per sample的中位数左右
--o-visualization alpha-rarefaction.qzv
--m-metadata-file samplemetadata.tsv \

#下载silva官方的silva-138-99-seqs.qza和silva-138-99-tax.qza
#训练分类器
#1，提取扩增区
qiime featue-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer CCTAYGGGRBGCASCAG \
--p-r-primer GGACTACMMGGGTATCTAAT \
--p-n-jobs 4 \
--o-reads ref-seqs.qza 

#2.只保留细菌的代表序列
qiime taxa filter-seqs \
--i-sequences ref-seqs.qza \
--i-taxonomy silva-138-99-tax.qza \
--p-include Bacteria \
--o-filtered-sequences ref-seqs-Bacteria.qza 

#3.仅保留细菌的注释信息，需要使用qiime2-2023.9-shotgun（服务器）
qiime rescript filter-taxa \
--i-taxonomy silva-138-99-tax.qza \
--m-ids-to-keep-file ref-seqs-Bacteria.qza \
--o-filtered-taxonomy ref-seqs-Bacteria-tax.qza 

#4.训练分类器
qiime featue-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs-Bacteria.qza \
--i-reference-taxonomy ref-seqs-Bacteria-tax.qza \
--o-classifier classifier-Bacteria.qza 

#物种注释（第一次）#需要1-2天
qiime featue-classifier classify-sklearn \
--i-reads rep-seqs.qza \
--i-classifier classifier-Bacteria.qza \
--o-classification taxonomy.qza 

#过滤特征表
#①过滤丰度<5的序列
qiime featue-table filter-featues \
--i-table table.qza \
--p-min-frequency 5 \
--o-filtered-table feature-frequency-filtered-table.qza 
#过滤线粒体和叶绿体序列
qiime taxa filter-table \
--i-table featue-frequency-filtered-table.qza \
--i-taxonomy taxonomy.qza \
--p-include p_ \
--p-exclude Mitochondria,chloroplast \
--o-filtered-table table-with-phyla-no-Mitochondria-no-chloroplast.qza
#更新特征序列
mkdir filtered_rep_seqs
qiime feature-table filter-seqs \
--i-data dada2/rep-seqs.qza \
--i-table filtered_table/table-with-phyla-no-Mitochondria-no-chloroplast.qza \
--o-filtered-data filtered_rep_seqs/final_rep_seqs.qza 
#导出特征序列

#物种注释（第二次）
qiime feature-classifier classify-sklearn \
--i-reads  final_rep_seqs.qza\
--i-classifier classifier-Bacteria.qza \
--o-classification final-taxonomy.qza 
#物种注释可视化
qiime metadata tabulate \
--m-input-file final-taxonomy.qza \
--o-visualization final-taxonomy.qzv 

#物种组成分析
绘制柱状图
qiime taxa barplot \
  --i-table table-with-phyla-no-Mitochondria-no-chloroplast.qza
  --i-taxonomy final-taxonomy.qza \
  --m-metadata-file samplemetadata.tsv \
  --o-visualization taxa-bar-plots.qzv  

#差异丰度分析(可作为组间比较的一种辅助分析的方法)
#进行ANCOM分析（ASV水平）
#对丰度表进行过滤，只比较盲肠/十二指肠样本组之间的差异，使其更符合ANCOM假设
qiime feature-table filter-samples \
--i-table table-with-phyla-no-Mitochondria-no-chloroplast.qza \
--m-metadata-file samplemetadata.txt \
--p-where "Posotion='duodenum'" \
--o-filtered-table duodenum-table.qza
#参数说明
## --i-table          输入特征丰度表
## --m-metadata-file  样本信息表,可使用元数据表
## --p-where          参数用双引号包围，[]中填写进行过滤要依据的列名，
##                    =号后面的''中填写该列名下进行过滤要保留的组名
## --o-filtered-table 输出过滤后的特征表

#1.特征丰度表添加伪计数，消除零值
qiime composition add-pseudocount \
--i-table duodenum-table.qza \
--o-composition-table comp-duodenum-table.qza 
#进行ANCOM分析（ASV水平）
qiime composition ancom \
  --i-table comp-duodenum-table.qza \
  --m-metadata-file samplemetadata.txt \
  --m-metadata-column group \
  --o-visualization ancom-duodenum-genetype.qzv


#2.按属水平合并统计
#对过滤后只保留肠道样本的特征丰度表进行汇总，获得属水平物种丰度表
qiime taxa collapse \
  --i-table duodenum-table.qza \
  --i-taxonomy final-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table duodenum-table-l6.qza
## 参数说明
##  --i-table           输入特征丰度表
##  --i-taxonomy        输入物种注释信息表
##  --p-level           要汇总的分类学水平
##  --o-collapsed-table 输出某分类学水平的物种丰度表 

# 丰度表添加伪计数，消除零值
qiime composition add-pseudocount \
  --i-table duodenum-table-l6.qza \
  --o-composition-table comp-duodenum-table-l6.qza 

#进行 ancom 分析
qiime composition ancom \
  --i-table comp-duodenum-table-l6.qza \
  --m-metadata-file samplemetadata.txt \
  --m-metadata-column group \
  --o-visualization l6-duodenum-ancom-genetype.qzv





#数据导出及R语言可视化




#qiime2导出物种丰度表
mkdir phyloseq
qiime tools export \
--input-path table.qza \#或者是过滤后的table-with-phyla-no-Mitochondria-no-chloroplast.qza
-output-path phyloseq  #biom格式，feature-table.biom

#qiime2+biom+qiime1获得16S物种丰度
#1.导出物种分类信息和置信度
qiime tools export \
--input-path final-taxonomy/final-taxonmony.qza --output-path taxa
#处理表头
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv
#导出otu(feature)表
qiime tools export \
   --input-path table-with-phyla-no-Mitochondria-no-chloroplast.qza \
   --output-path table_exported
#进入linux
biom add-metadata \
-i feature-table.biom \
-o featue-table_w_tax.biom \
--observation-metadata-fp taxa/taxonomy.tsv
--sc-separated taxonomy
#进入ubuntu,激活qiime1
#结果按门、纲、目、科、属五个级别进行分类汇总，对应结果的L2-L6
summarize_taxa.py -i feature-table_w_tax.biom -o result/sum_taxa # summary each level percentage，相对丰度表
#结果导入excel表格中手动添加metadata分组信息,格式如下：
#bodysite                                mucosal         mucosal         mucosal         mucosal         mucosal         non_mucosal     non_mucosal     non_mucosal     non_mucosal     non_mucosal
#subsite                                 oral            gut             oral            oral            gut             skin            nasal           skin            ear             nasal
#id                                      1023            1023            1672            1876            1672            159005010       1023            1023            1023            1672
#Bacteria                                0.99999         0.99999         0.999993        0.999989        0.999997        0.999927        0.999977        0.999987        0.999997        0.999993
#Bacteria|Actinobacteria                 0.311037        0.000864363     0.00446132      0.0312045       0.000773642     0.359354        0.761108        0.603002        0.95913         0.753688
#Bacteria|Bacteroidetes                  0.0689602       0.804293        0.00983343      0.0303561       0.859838        0.0195298       0.0212741       0.145729        0.0115617       0.0114511
#Bacteria|Firmicutes                     0.494223        0.173411        0.715345        0.813046        0.124552        0.177961        0.189178        0.188964        0.0226835       0.192665
#Bacteria|Proteobacteria                 0.0914284       0.0180378       0.265664        0.109549        0.00941215      0.430869        0.0225884       0.0532684       0.00512034      0.0365453
#Bacteria|Firmicutes|Clostridia          0.090041        0.170246        0.00483188      0.0465328       0.122702        0.0402301       0.0460614       0.135201        0.0115835       0.0537381

#数据格式转换
#激活lefse
#-c指定分组行；-s指定亚组行，若没有可以不指定；-u指定样本编号；-o指定归一化后范围，主要针对宏基因组数据，目的是对相对丰度进行放大
lefse_format_input.py hmp_aerobiosis_small.txt hmp_aerobiosis_small.in -c 1 -s 2 -u 2 -o 1000000
#lefse分析
lefse_run.py hmp_aerobiosis_small.in hmp_aerobiosis_small.res
参数设置：-a指定组间比较检验水准阈值，-w指定成组比较检验水准阈值，-l指定lda score阈值。除了可以选择lda，还可以选择svm进行分析。
#绘制lefse结果图
lefse_plot_res.py hmp_aerobiosis_small.res hmp_aerobiosis_small.png --dpi 500 



ubuntu可能会出现无法显示共享文件夹的bug
在用户的主目录下使用，sudo vmhgfs-fuse .host:/ /mnt/hgfs/ -o allow_other -o uid=1000

#如何得到没有零值的相对丰度表
将过滤后的table.qza添加零值即可（跟ANCOM去除零值相同）

picrust2功能预测
### 自己手动运行#链接文件过来，qiime2分析的table和序列ln -s ../4.export_data/feature_table_tax.biom .ln -s ../4.export_data/dna-sequences.fasta .  ##两个文件：dna-sequences.fasta featue-table.biom
#PICRUSt2提供流程一键完成f分析picrust2_pipeline.py -s dna-sequences.fasta -i featue-table.biom -o picrust2_out_pipeline -p 20
#后续分析结果整理metagenome_pipeline.py -i featue-table.biom -m picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz -f picrust2_out_pipeline/EC_predicted.tsv.gz \                       
-o EC_metagenome_out 
--strat_outadd_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \                   
 -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gzconvert_table.py EC_metagenome_out/pred_metagenome_contrib.tsv.gz \                 
 -c contrib_to_legacy \                 -o EC_metagenome_out/pred_metagenome_contrib.legacy.tsv.gz

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \                    
-o pathways_out -p 10add_descriptions.py 

pathway_pipeline.py -i pathways_out/path_abun_unstrat.tsv.gz 
-m METACYC \                   
-o pathways_out/path_abun_unstrat_descrip.tsv.gz

#结果注释
在第一列后面添加一列注释，结果表也方便在STAMP中进行差异比较。
添加EC的注释
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
KO添加注释
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
pathway添加注释，基于MetaCyc数据库
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o pathways_out/path_abun_unstrat_descrip.tsv.gz





##一个基于功能预测数据下游分析的R包ggpicrust2
安装R包
BiocManager::install("lefser")
devtools::install_github('cafferychen777/ggpicrust2')
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
