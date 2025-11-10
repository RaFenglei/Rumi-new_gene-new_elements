# 使用ANNEVO来注释基因组
[ANNEVO](https://github.com/xjtu-omics/ANNEVO)是一种基于深度学习的从头开始的基因注释方法，用于理解基因组功能。ANNEVO 能够直接从基因组中模拟不同物种的远端序列信息和联合进化关系。


## Step 1：安装
比较推荐两种方式
1.使用git clone安装然后编译
```sh
git clone https://github.com/xjtu-omics/ANNEVO.git
cd ANNEVO
```
2.直接使用conda安装
```sh
conda create -n ANNEVO_v2 python=3.10
```
**需要注意的是，我们服务器暂时没有配置GPU，所以安装好直接使用即可，无需再安装CUDA等配套软件**

## Step 2：使用
**由于我们没有GPU，所以只能使用One-step Execution直接生成gff，经过测试，一个2.5G左右的哺乳动物基因组使用哺乳动物模型来执行的话，CPU大约需要1天，GPU大约需要5分钟**
```sh
echo Start Time is `date`
cd /data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/02.ANNEVO/01.deer_test

#set -e
#source ~/bin/cactus-bin-v1.2.3/venv/bin/activate

/public/home/liyongxin/lilisen/miniconda3/envs/ANNEVO/bin/python /public/home/liyongxin/lilisen/soft/ANNEVO/ANNEVO-main/annotation.py --genome mhl.hifiasm.hic.hap1.fasta --lineage Mammalia --output mhl.annevo.gff --threads 128

echo End Time is `date`
```
**但是如果你是有GPU的，还是建议使用下面的Step-by-step Execution生成**
1：预测每个核苷酸的三种信息类型（建议在 GPU 资源丰富的环境中执行）
2：将 3 类信息解码为生物学上有效的基因结构（建议在 CPU 资源丰富的环境中进行）
```sh
# Nucleotide prediction
python prediction.py --genome path_to_genome --model_path path_to_model --model_prediction_path path_to_save_predction

# Gene structure decoding
python decoding.py --genome path_to_genome --model_prediction_path path_to_save_predction --output path_to_gff --threads 48 
```

## Step 3：重训练模型
官方只提供了Mammalia等比较常见的谱系模型，如果需要在特定进化枝上合并其他物种或重新训练 ANNEVO （比如利用注释最好的牛羊鹿训练反刍谱系模型），可以按照以下脚本进行作：
```sh
train_species_list="The species list used for training model"
val_species_list="The species list used for validating model"
h5_data_path="The path to store h5 file" 
mkdir -p tmp

# The file must be cleared before each run.
rm -f ${h5_data_path}/train.h5 ${h5_data_path}/train_with_intergenic.h5
rm -f ${h5_data_path}/val.h5 ${h5_data_path}/val_with_intergenic.h5

for species_name in "${train_species_list[@]}"; do
    path_to_genome="The path to species genome"
    path_to_annotation="The path to species annotation"
    # Filter out duplicated gene IDs and other issues that may cause parsing errors in the Biopython package
    python src/filter_wrong_record.py --input_file ${path_to_annotation} --output_file "tmp/tmp_${species_name}.gff"
    # Convert the genome sequence and annotation into H5 data for model training.
    python generate_datasets.py --genome ${path_to_genome} --annotation "tmp/tmp_${species_name}.gff" --output_file "${h5_data_path}/train" --threads 64
    rm -f "tmp/tmp_${species_name}.gff"
done
for species_name in "${val_species_list[@]}"; do
    path_to_genome="The path to species genome"
    path_to_annotation="The path to species annotation"
    python src/filter_wrong_record.py --input_file ${path_to_annotation} --output_file "tmp/tmp_${species_name}.gff"
    python generate_datasets.py --genome ${path_to_genome} --annotation "tmp/tmp_${species_name}.gff" --output_file "${h5_data_path}/val" --threads 64
    rm -f "tmp/tmp_${species_name}.gff"
done

# Train the deep learning model
python model_train.py --h5_path ${h5_data_path} --model_save_path path_to_new_model.pt
```

**当然，如果只是发现与你的物种密切相关的物种的基因组有限或者不可使用时间，可以选择 ANNEVO 的五种主要训练模型之一作为微调的起点**
```sh
# Filter out duplicated gene IDs and other issues that may cause parsing errors in the Biopython package
fine_tune_species_list="The species list used for fine tuning model"
h5_data_path="The path to store h5 file"
mkdir -p tmp

# The file must be cleared before each run.
rm -f ${h5_data_path}/fine_tune.h5 ${h5_data_path}/fine_tune_with_intergenic.h5

for species_name in "${fine_tune_species_list[@]}"; do
    path_to_genome="The path to species genome"
    path_to_annotation="The path to species annotation"
    python src/filter_wrong_record.py --input_file ${path_to_annotation} --output_file "tmp/tmp_${species_name}.gff"
    python generate_datasets.py --genome ${path_to_genome} --annotation "tmp/tmp_${species_name}.gff" --output_file "${h5_data_path}/fine_tune" --threads 64
    rm -f "tmp/tmp_${species_name}.gff"
done

# Fine tuning deep learning model
python fine_tune.py --model_path path_to_existing_model.pt --model_save_path path_to_new_model.pt --h5_path ${h5_data_path}
```

## Step 4：BUSCO验证结果质量





