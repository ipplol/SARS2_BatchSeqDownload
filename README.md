# SARS2_BatchSeqDownload
This script is designed to automatically download all sequence mutations and corresponding metadata from the [RCoV19](https://ngdc.cncb.ac.cn/ncov/?lang=en) database which includes data from GISAID, GenBank, and Other sources.

# Data format
No fasta file is available. Instead, it directly downloads the vcf file and converts it to a lighter and more readable format called the **sample_mutlist.tsv**. It is a table-separated file containing two columns, the ACCESSION ID and the list of mutations which is further separated by SPACE.
*"OU161829.1	7771T/C 26144G/T 17247T/C 3851G/A 1515A/G 14805C/T 12565G/T"*

This file will be further combined with metadata downloaded like the Lineage, Collection data, and the Collection location. Generating a file like this:

*"EPI_ISL_7489404	210G/T 241C/T 3037C/T 4181G/T 5693C/T 6402C/T 7124C/T 7482C/T 7851C/T 8835T/C 8986C/T 9053G/T 10029C/T 10969C/T 11201A/G 11332A/G 12513C/T 14408C/T 15451G/A 15521T/A 16466C/T 19220C/T 21618C/G 21846C/T 21987G/A 22028GAGTTCA/G 22082C/T 22858C/T 22917T/G 22995C/A 23403A/G 23604C/G 24410G/A 25469C/T 26767T/C 27638T/C 27752C/T 27874C/T 28247AGATTTC/A 28270TA/T 28461A/G 28881G/T 28916G/T 29402G/T 29742G/T	United Kingdom	2021-12-03	AY.4"*

Finally, a script will annotate the mutation in the RBD region, generating another file like:

*"EPI_ISL_7489404	L452R T478K	United Kingdom	2021-12-03	AY.4"*

Due to the massive number of sequences, this script can require up to 600GB of **empty storage space** (as of 2023.02.14). 

# Command
The command is very simple.

Once you have downloaded the whole files from GitHub.
Go to the script folder.
**"cd ./BSD/Script"**
And bash the shell.
**"bash DataPrepare.sh"**

It will automatically start to compile the c++ scripts and then start downloading and data processing. **All files downloaded will be stored in the ./BSD/Data folder.**

Please notice that the whole process can take up to 12 hours to complete depending on the performance of your computer.
