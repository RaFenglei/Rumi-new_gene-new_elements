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
cat > ./01.Cervidae_deer_ref/01.hal_to_maf.sh << 'EOF'
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
EOF
```

投递运行完子脚本后，我们就可以在./01.Cervidae_deer_ref/deerref找到


## Step 3：清理生成的maf文件
在cactues-hal2maf得到的结果虽然我们去除了其中的Duplicate，但是仍然会存在一比多或者多比一的情况，所以我们需要对生成的maf文件进行清理才能进行后续的操作

```sh
cat > ./01.Cervidae_deer_ref/deerref/01_1.mafclean.sh << 'EOF'
#usage : sh mafclean.sh [maf] [speicelist]
maf=$1
sp_l1=$(cat $2|xargs|sed 's/ /,/g')
sp_l2=$(cat $2|xargs)
#maf_parse -o MAF -n -O "$sp_l2" "$maf"
#printf "maf_order  <(maf_parse -o MAF -n -O ${sp_l1} ${maf}) ${sp_l2} > $maf.clean.maf"
maf_order  <(maf_parse -o MAF -n -O ${sp_l1} ${maf}) ${sp_l2} > $maf.clean.maf
EOF
```

当然，因为我们一般是按照染色体分割的，实在是太多了，所以可以直接下面一行命令把所有的maf全部处理完

```sh
cat > ./01.Cervidae_deer_ref/deerref/01_2.mafclean.sh << 'EOF'
source /public/home/wangwen_lab/lizihe/.bashrc
find /data02/zhangfenglei/project/09.new_gene_elements/03.new_peak/01.hal2maf/01.verson1/01.Cervidae_deer_ref/deerref -name "*.maf" -exec bash -c 'for f; do bash 02_1.mafclean.sh "$(basename "$f")" species_list; done' _ {} +
EOF
```


## Step 4：排序生成的maf文件
虽然我们已经去除了maf中的Duplicate并且剔除了其中的异常比对情况，但是不管是计算HCE还是Syntenic Region，所以还需要对生成的maf进行排序
我们运行下面的02.mafsort.sh即可
```sh
cat > ./01.Cervidae_deer_ref/deerref/02.mafsort.sh << 'EOF'
source /public/home/wangwen_lab/lizihe/.bashrc
find /data02/zhangfenglei/project/09.new_gene_elements/03.new_peak/01.hal2maf/01.verson1/01.Cervidae_deer_ref/deerref -name "*.clean.maf" -exec bash -c 'for f; do /public/home/wangwen_lab/zhoubotong/soft/last-1061/scripts/maf-sort "$f" > "${f%.clean.maf}.clean.maf.sort.maf"; done' _ {} +
EOF
```

你需要注意的是，不管是maf的clean还是sort，都需要一个species_list，这里面每一行一个你需要的参考物种，并且注意顺序就是你最终得到的.clean.maf.sort.maf中的物种排序
```species_list
cat > species_list << 'EOF'
mhl
ML
PL
TUNLU
TL
ECEP
reev_muntjac
PZ
mule_deer
XL
goat
sheep
cattle
muskdeer
okapi
giraffe
pronghorn
xilu
hippo
gray_whale
killer_whale
pig
camel
tapir
human
mouse
EOF
``` 



## Reference:

[Cactus论文](https://www.nature.com/articles/s41586-020-2871-y)

[实操：蝴蝶论文](https://www.nature.com/articles/s41467-024-50529-0)
