# Panda姐的宏基因组入门-工作环境搭建（1）

### 1 Minconda软件管理中心

参照https://github.com/shenmengyuan/RNA_seq_Biotrainee中的1.2 安装minconda

```
# 下载安装包（Linux版本）
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# 也可以在清华镜像下载：https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/ 
# wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 直接安装,一路yes即可
bash Miniconda3-latest-Linux-x86_64.sh

# 添加生信软件包下载频道
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# 下面是清华的频道镜像：https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
# 据说有时候不稳定，大家可以尝试下：
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
# conda config --set show_channel_urls yes
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/menpo/
```

### 2 anvio工作环境

```shell
conda create -n anvio4 -c bioconda -c conda-forge python=3 anvio=4
# 大概需要下载安装714.0 MB的软件
source activate anvio4 # 进入工作环境
source deactivate # 离开工作环境
```

需要安装的软件有：

- `bwa-0.7.17`（用于比对，也可以自己安装`bowtie2`）
- `diamond-0.9.19`（替代`blastp`,`blastx`）
- `muscle-3.8.1551`（用于多重序列比对，常见于使用蛋白质进化树的构建时的比对）
- `hmmer-3.1b2`（基于隐马尔科夫模型的比对软件）
- `prodigal-2.6.3`（基因预测）
- `blast-2.2.31 `（用于比对）
- `samtools-1.4.1`（操作`bam`,`sam`）
- `centrifuge-1.0.3`（用于序列物种鉴定）
- `anvio-4.0.0`（contig可视化）
- `checkm`（用于分箱或者微生物基因组拼接质控）
- `metabat`（用于contig分箱）
- `eggnog`（用于`go`,`cog`注释）
- `megan6`（用于功能注释、分类注释，同为`diamond`开发的团队开发）
- `kaiju`（基于kmer的快速物种分类）

前面九个软件在`anvio`的`conda`工作环境里有了，下面是手动安装的软件：

```shell
## 安装checkm
# 官网：http://ecogenomics.github.io/CheckM/
conda create -n checkm_python2.7 python=2.7
source activate checkm_python2.7
pip install checkm-genome

## 安装metabet
# 官网：https://bitbucket.org/berkeleylab/metabat
# 下载编译好的包，下载解压即可用；
wget https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v2.12.1.tar.gz

## 安装eggnog
# 官网：http://eggnogdb.embl.de/
wget  https://github.com/jhcepas/eggnog-mapper/archive/1.0.1.tar.gz
tar xvzf 1.0.1.tar.gz

## 安装megan
# 官网：http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html
axel http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_11_0.sh

## 安装kaiju
# 官网：http://kaiju.binf.ku.dk/
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make


```

### 3 数据库

```shell
## centrifuge index数据库
# NCBI nucleotide non-redundant sequences 
axel ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/nt_2018_2_12.tar.gz
# Bacteria, Archaea (compressed)
axel ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed.tar.gz

## eggnog对应的细菌数据库下载
download_eggnog_data.py bact

## silva 原核和真核微生物的小亚基rRNA基因序列（简称SSU，即16S和18SrRNA)
axel https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz

## nr非冗余数据库
# 结合diamond进行nr库比对
axel ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
axel ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5
# https://github.com/bbuchfink/diamond
# diamond makedb --in nr.faa -d nr

## MEGAN注释文件
axel http://ab.inf.uni-tuebingen.de/data/software/megan6/download/prot_acc2tax-Mar2018X1.abin.zip
axel http://ab.inf.uni-tuebingen.de/data/software/megan6/download/SSURef_NR99_128_tax_silva_to_NCBI_synonyms.map.gz

## kaiju物种注释文件
# Representative genomes from proGenomes
makeDB.sh -p -v
# Non-redundant protein database nr
makeDB.sh -n 
```



### 4 参考资料

[微生物组软件导读](http://mp.weixin.qq.com/s/cVVp2YEVKwg1bcTgCEj8sw)

[CheckM：微生物基因组组装/metagenome基因组重构完整度和杂合度评估](https://mp.weixin.qq.com/s/3tfsN32lIuF08YJ9_Wclbw)

[Centrifuge：又是一款快速有效的 metagenome 序列分类的软件](http://mp.weixin.qq.com/s/8VbfaWj2pcHgOx91uQv9Ww)

[MEGAN6的安装和使用笔记](http://mp.weixin.qq.com/s/SewGe-hwoF1Ft10INlR7mw)

[diamond安装及使用说明阅读笔记](http://mp.weixin.qq.com/s/LTLX1a3YxLhGgFhdRTSyTA)

