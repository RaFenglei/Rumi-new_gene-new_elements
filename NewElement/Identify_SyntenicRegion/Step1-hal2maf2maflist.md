# 使用Cactus得到多物种比对的maf文件，然后转化为maflist文件作为输入来鉴定保守的插入和缺失
[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/tree/master)是一个无参考的全基因组比对程序，并且也可以用来构建Pan-genome
**本文以新基因新元件课题为例，展示了如何从Cactus得到的hal文件中一步一步得到maf文件和maflist文件**

参考教程一：从hal中提取maf文件，https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/maf/cactus_hal2maf.py

参考教程二：使用本方法发表的文献，https://www.nature.com/articles/s41467-024-50529-0

参考路径：
hal2maf参考路径/data02/zhangfenglei/project/09.new_gene_elements/03.new_peak/01.hal2maf/01.verson1

## Step 1：安装
Cactus安装请参考github，本文直接使用了Bao Wang(Postdoc, Northwestern Polytechnical University)和Zecheng Du(PhD Candidate, Northwestern Polytechnical University)路径下的Cactus
提取现有HAL文件的根节点


## Step 2：从HAL文件中提取对应的maf文件
我们首先使用halStats查看hal中的根节点，发现是Anc00
```sh

mkdir -p qsub_shell
mkdir -p deerref
mkdir -p tmp

declare -a chromosomes=("chr1a" "chr2a" "chr3a" "chr4a" "chr5a" "chr6a" "chr7a" "chr8a" "chr9a" "chr10a" "chr11a" "chr12a" "chr13a" "chr14a" "chr15a" "chr16a" "chr17a" "chr18a" "chr19a" "chr20a" "chr21a" "chr22a" "chr23a" "chr24a" "chr25a" "chr26a" "chr27a" "chr28a" "chr29a" "chr30a" "chr31a" "chr32a" "chrX" "chrY")

for i in "${!chromosomes[@]}"; do
    chrom_name="${chromosomes[i]}"
    echo "source /public/home/wangwen_lab/zhouxiaofang/.bashrc
cd /data02/zhangfenglei/project/09.new_gene_elements/03.new_peak/01.hal2maf/01.verson1/01.Cervidae_deer_ref
source activate python3
source ~zhoujiong/software/cactus-bin-v2.8.2/venv-cactus-v2.8.2/bin/activate
cd deerref
cactus-hal2maf ./js$(($i+1)) /data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactus/04.add_xilu_new/01.existing_tree.hal 26sp.${chrom_name}.deer_ref.hal2maf.maf.gz --chunkSize 500000 --refGenome mhl --refSequence $chrom_name --noAncestors --batchCores 32 --filterGapCausingDupes --targetGenomes mhl,ML,PL,TUNLU,TL,ECEP,reev_muntjac,PZ,mule_deer,XL,goat,sheep,cattle,muskdeer,okapi,giraffe,pronghorn,xilu,hippo,gray_whale,killer_whale,pig,camel,tapir,human,mouse --workDir /data02/zhangfenglei/project/09.new_gene_elements/03.new_peak/01.hal2maf/01.verson1/01.Cervidae_deer_ref/tmp
zcat 26sp.${chrom_name}.deer_ref.hal2maf.maf.gz | mafDuplicateFilter -k -m - > 26sp.${chrom_name}.deer_ref.hal2maf.single.maf" > qsub_shell/mhl_$(($i+1)).sh
done
```

然后
我们选择的物种是mhl,ECEP,PZ,XL,sheep,cattle,muskdeer,giraffe,pronghorn,xilu作为反刍动物每个科（其中鹿科是每个族、牛科是每个亚科选了一个代表物种）来得到相对于其他所有物种的(mhl,ML,PL,TL,TUNLU,ECEP,reev_muntjac,HJ,PZ,mule_deer,XL,goat,sheep,cattle,muskdeer,okapi,giraffe,pronghorn,xilu,hippo,gray_whale,killer_whale,pig,camel,tapir,human)chain文件


## Step 3：创建新的系统发育树定义
在现有树(Anc00:0.0162236,tapir:0.0162236)Anc00;的基础上，添加human分支
需要注意的是，0.0162236只是预估值，填入进去cactus还会重新计算
```sh
cat > 02.cactus_add_human.in << 'EOF'
(Anc00:0.0162236,human:0.0162236)AncHuman;
Anc00 /data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactuss/02.add_human/01.existing_tree.fa
human /data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactuss/02.add_human/00.human.fa
EOF
```


## Step 4：创建Cactus作业脚本
```sh
cat > 02.cactus_add_human.sh << 'EOF'
#!/usr/bin/bash
#PBS -V
set -exo

source /public/home/wangwen_lab/lizihe/.bashrc
conda activate /public/home/wangwen_lab/lizihe/soft/anaconda3/envs/cactus
export PATH=/public/home/wangwen_lab/lizihe/soft/cactus-bin-v2.9.4/bin:$PATH

cd /data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactuss/02.add_human
mkdir -p jobstore_human jobstore_human/logs
cp 02.cactus_add_human.in jobstore_human/

## Preprocessor
cactus-preprocess ./jobstore_human/0 02.cactus_add_human.in jobstore_human/02.cactus_add_human.in \
    --maxCores 100 --inputNames Anc00 human \
    --logFile jobstore_human/logs/preprocess-AncHuman.log

## Alignment
### Blast
cactus-blast ./jobstore_human/1 jobstore_human/02.cactus_add_human.in jobstore_human/AncHuman.paf \
    --root AncHuman --maxCores 100 \
    --logFile jobstore_human/logs/blast-AncHuman.log

### Align
cactus-align ./jobstore_human/2 jobstore_human/02.cactus_add_human.in jobstore_human/AncHuman.paf jobstore_human/out_with_human.hal \
    --root AncHuman --maxCores 100 \
    --logFile jobstore_human/logs/align-AncHuman.log

## 提取最终的human祖先序列（可选）
hal2fasta jobstore_human/out_with_human.hal AncHuman --hdf5InMemory > jobstore_human/AncHuman.fa

## 验证结果
halValidate jobstore_human/out_with_human.hal --hdf5InMemory

EOF
```

**最终输出结构是**
```txt
/data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactuss/02.add_human/
├── jobstore_human/                    # 主要输出目录
│   ├── 02.cactus_add_human.in         # 复制的输入文件
│   ├── AncHuman.paf                   # Blast结果文件
│   ├── out_with_human.hal             # **最终输出的HAL文件**
│   ├── AncHuman.fa                    # 祖先序列文件（可选）
│   └── logs/                          # 日志目录
│       ├── preprocess-AncHuman.log
│       ├── blast-AncHuman.log
│       └── align-AncHuman.log
├── existing_tree.fa                   # 现有树的根节点序列
└── 02.cactus_add_human.in            # 原始输入文件
```
jobstore_human/out_with_human.hal​​ - ​​最重要的文件​​，包含human的完整多序列比对结果

​​jobstore_human/AncHuman.paf​​ - Blast比对结果

​​jobstore_human/AncHuman.fa​​ - 新根节点AncHuman的祖先序列

同时你需要注意的是
（1）如果是在计算节点直接运行，请先加载环境
``` sh
source /public/home/wangwen_lab/lizihe/.bashrc
conda activate /public/home/wangwen_lab/lizihe/soft/anaconda3/envs/cactus
export PATH=/public/home/wangwen_lab/lizihe/soft/cactus-bin-v2.9.4/bin:$PATH

bash 02.cactus_add_human.sh > 02.cactus_add_human.log 2>&1
```
（2）如果是投递到qsub运行，直接投递就好
``` sh
#source /public/home/wangwen_lab/lizihe/.bashrc
#conda activate /public/home/wangwen_lab/lizihe/soft/anaconda3/envs/cactus
#export PATH=/public/home/wangwen_lab/lizihe/soft/cactus-bin-v2.9.4/bin:$PATH
qsub -q high1 -l nodes=1:ppn=128 02.cactus_add_human.sh
``` 



## Reference:

[Cactus论文](https://www.nature.com/articles/s41586-020-2871-y)

[Hal论文](https://academic.oup.com/bioinformatics/article/29/10/1341/256598?login=false)
