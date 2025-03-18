# 这是一个自用生信分析脚本
## 如何运行
1. 将文件下载，放在同一个目录
2. 确保docker已经开启并已登录
3. 运行build.bat
4. 将需要分析的文件放在input文件夹，默认分析.fastq.gz文件，如果要分析.fasta文件，请打开pipeline.sh，修改`your_input_file_type_is`变量的字符`"fastq"`为`"fasta"`
5. 将参考基因组文件更名为reference.fasta并放在reference文件夹
6. 运行run.bat
## 做了什么
基于docker的ubuntu镜像，利用fastp、SPAdes、quast、snippy、gubbins工具完成质量控制、序列拼接、拼接质量评估、突变位点分析和构建进化树。
## 怎么看结果
所有结果文件都在output文件夹，不同插件的输出文件均在其同名文件夹中，SPAdes拼接出来的序列会在fasta文件夹中。