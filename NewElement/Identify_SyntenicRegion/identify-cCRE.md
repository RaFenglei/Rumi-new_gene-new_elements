# 使用Cactus得到两两物种比对的chain文件，然后转化为Reciprocal Best的chain作为TOGA的输入来鉴定新基因
[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/tree/master)是一个无参考的全基因组比对程序，并且也可以用来构建Pan-genome
**本文以新基因新元件课题为例，展示了如何从Cactus得到的hal文件中一步一步得到两两物种比对的Reciprocal Best Chian文件**

参考教程一：从hal中提取chain文件，https://ucsc-ci.com/comparativegenomicstoolkit/cactus/-/blob/167179fcff6e09d7864170013c15029da04d8841/doc/progressive.md

参考教程二：处理chain文件得到Reciprocal Best Chian文件，http://genomewiki.ucsc.edu/index.php/HowTo:_Syntenic_Net_or_Reciprocal_Best

参考路径：/data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactus/03.hal2chain

## Step 1：安装
Cactus安装请参考github，本文直接使用了Bao Wang(Postdoc, Northwestern Polytechnical University)和Zecheng Du(PhD Candidate, Northwestern Polytechnical University)路径下的Cactus
提取现有HAL文件的根节点


## Step 2：提取现有HAL文件的所以genome
我们首先使用halStats查看hal中的根节点，发现是Anc00
```txt
$halStats out_with_human.hal

hal v2.2
(((camel:0.0244494,(pig:0.0269299,((xilu:0.0318519,((pronghorn:0.0119009,(giraffe:0.00409218,okapi:0.00623621)Anc10:0.00487476)Anc08:0.000706151,((muskdeer:0.00961721,(cattle:0.0071056,(sheep:0.00238328,goat:0.0023734)Anc16:0.0058297)Anc13:0.00197833)Anc11:0.000900202,(((reev_muntjac:0.0041239,ECEP:0.00390819)Anc17:0.000983137,(TUNLU:0.00275209,(TL:0.00249469,((mhl:0.00071661,ML:0.000762197)Anc23:0.000915221,PL:0.00192979)Anc22:0.000477065)Anc21:0.000227643)Anc18:0.00125769)Anc14:0.000957613,((HJ:0.00410087,PZ:0.00410087)Anc19:0.00263453,(mule_deer:0.00285365,XL:0.00281675)Anc20:0.00210855)Anc15:0.000927252)Anc12:0.0056893)Anc09:0.000885046)Anc06:0.0111833)Anc04:0.0106196,(hippo:0.0192975,(gray_whale:0.00756346,killer_whale:0.00970889)Anc07:0.0108208)Anc05:0.00218731)Anc03:0.00439133)Anc02:0.00331519)Anc01:0.0162236,tapir:0.0162236)Anc00:0.0162236,human:0.0162236)AncHuman;
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
