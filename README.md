# DOWAS manual
## update date: 2018/10/17
## version 1.0.3
## Getting started
```
cd /your/path
tar czvf DOWASv1.0.3.tar.gz 
cd DOWAS v1.0.3
perl run.pl example.xml thread_num
```
## Introduction
DOWAS is a perl script for find potential off-target sites from reference and germline sequencing data in human. It consists of three pl scripts, run.pl will invoke other 2 scripts. All information you need to input should be written in a xml file like example.xml and run ‘perl run.pl example.xml 12’ in your computer. in that case, DOWAS will generate a folder in current path with same name with xml file (in that case it names ‘example’) and put all results in this folder and this software will use 12 threads. 
You need to install some software and write them in the xml file before using DOWAS, those softwares with recommend version are list at the end of this manual.
Xml file contain a ‘rawdata_list’ subtag, it at least contains a father tag and a mother tag which contain raw sequencing data from parents in germline.
Final result will be saved in /your/path/example/$offspring/final_result.txt, you will get multiple offspring folder and multiple final_result.txt if you have more than one offspring.
## Availability
DOWAS is released under MIT license. The latest version of DOWAS can be downloaded at www.sustc-genome.edu.cn. 
## Seeking help
If you have questions about DOWAS, you may send the questions to chenyr@mail.sustc.edu.cn or chenkj@mail.sustc.edu.cn or 11749245@mail.sustc.edu.cn or caow@mail@sustc.edu.cn . You may also ask questions in forums such as BioStar and SEQanswers.
## Author
DOWAS is written by Chen Kaijing, Chen Yangran, ChenRui and Cao Wei.

## Reference
•	Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] 
•	Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
•	Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
•	Poplin R, Ruano-Rubio V, DePristo M A, et al. Scaling accurate genetic variant discovery to tens of thousands of samples[J]. bioRxiv, 2017: 201178.
## Software and database dependencies
fastp 0.19.4
samtools-1.9
jre-1.8.0_181
GATK-3.8
GATK-4.0.4.0
bwa-0.7.17
dbsnp database
OneKG_indel database
Mills_indel database
Human genome sequence reference
## update info
### DOWAS v1.0.1:
We remove the dependency software in DOWAS package to reduce the size of DOWAS and avoid copyright problem, now you need to download those specified version softwares and write address of them in configuration xml file. 
### DOWAS v1.0.2: 
We changed sequence extract function in our code to reduce about one day’s DOWAS running time.
### DOWAS v1.0.3:
We add a thread number parameter to let user choose a suitable thread number to run DOWAS. And we use pipeline in mapping period to reduce IO spending.
