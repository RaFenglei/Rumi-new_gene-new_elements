# 使用Cactus在最外枝添加基因组
[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/tree/master)是一个无参考的全基因组比对程序，并且也可以用来构建Pan-genome
**本文以新基因新元件课题为例，展示了如何在Cactus现有系统发育树最外层添加一个human分支**

参考路径：/data02/zhangfenglei/project/09.new_gene_elements/02.new_gene/01.cactuss/02.add_human

## Step 1：安装
Cactus安装请参考github，本文直接使用了Bao Wang(Postdoc, Northwestern Polytechnical University)和Zecheng Du(PhD Candidate, Northwestern Polytechnical University)路径下的Cactus
提取现有HAL文件的根节点


## Step 2：提取现有HAL文件的根节点
我们首先使用halStats查看hal中的根节点，发现是Anc00
```txt
$halStats ../01.genome/out.hal

hal v2.2
((camel:0.0244494,(pig:0.0269299,((xilu:0.0318519,((pronghorn:0.0119009,(giraffe:0.00409218,okapi:0.00623621)Anc10:0.00487476)Anc08:0.000706151,((muskdeer:0.00961721,(cattle:0.0071056,(sheep:0.00238328,goat:0.0023734)Anc16:0.0058297)Anc13:0.00197833)Anc11:0.000900202,(((reev_muntjac:0.0041239,ECEP:0.00390819)Anc17:0.000983137,(TUNLU:0.00275209,(TL:0.00249469,((mhl:0.00071661,ML:0.000762197)Anc23:0.000915221,PL:0.00192979)Anc22:0.000477065)Anc21:0.000227643)Anc18:0.00125769)Anc14:0.000957613,((HJ:0.00410087,PZ:0.00410087)Anc19:0.00263453,(mule_deer:0.00285365,XL:0.00281675)Anc20:0.00210855)Anc15:0.000927252)Anc12:0.0056893)Anc09:0.000885046)Anc06:0.0111833)Anc04:0.0106196,(hippo:0.0192975,(gray_whale:0.00756346,killer_whale:0.00970889)Anc07:0.0108208)Anc05:0.00218731)Anc03:0.00439133)Anc02:0.00331519)Anc01:0.0162236,tapir:0.0162236)Anc00;
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


## Reference:

[Cactus论文](https://www.nature.com/articles/s41586-020-2871-y)

[Hal论文](https://academic.oup.com/bioinformatics/article/29/10/1341/256598?login=false)
