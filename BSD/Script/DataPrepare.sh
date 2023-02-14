#script for download SARS-CoV-2 vcf and metadata file from the RCoV19 database
#and transfer them to a human readable format (tsv)
#!/bin/bash

echo "Compiling [Some warnings may show up which are normal, please relax.]:"
echo "Compiling [Some warnings may show up which are normal, please relax.]:"
echo "Compiling [Some warnings may show up which are normal, please relax.]:"
echo "Compiling [Some warnings may show up which are normal, please relax.]:"
echo "Compiling [Some warnings may show up which are normal, please relax.]:\n\n\n"

thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#compiling c++
cd `dirname "thisfold"`
g++ $thisfold/SeqFindMeta/main.cpp -o $thisfold/SeqFindMeta.out -std=c++11 -lpthread
g++ $thisfold/Mutresult2Sample_mutlistLINUX/Mutresult2Sample_mutlistLINUX/main.cpp -o $thisfold/Mutresult2Sample_mutlistLINUX.out -std=c++11 -lpthread
g++ $thisfold/Vcf2mut5col/Vcf2mut5col/main.cpp -o $thisfold/Vcf2mut5col.out -std=c++11 -lpthread
g++ $thisfold/Mut2Haplo/main.cpp -o $thisfold/Mut2Haplo.out -std=c++11 -lpthread

#Data update
mainDir1=`dirname $thisfold`
mainDir=`dirname $mainDir1`
cd "$mainDir/BSD/Data"
#wget download
wget ftp://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz -O 2019-nCoV_total.vcf.gz
wget https://ngdc.cncb.ac.cn/ncov/genome/export/meta -O meta
mv meta meta.zip
unzip meta.zip
gunzip 2019-nCoV_total.vcf.gz

#data process
echo "Data processing. It may take a while."
$mainDir/BSD/Script/Vcf2mut5col.out $mainDir/BSD/Data/2019-nCoV_total.vcf $mainDir/BSD/Data
echo "Done.\nStart Convert Mutresult to Sample_mutlist"
$mainDir/BSD/Script/Mutresult2Sample_mutlistLINUX.out $mainDir/BSD/Data
echo "Done.\nStart Combine Sample_mutlist and Metadata"
$mainDir/BSD/Script/SeqFindMeta.out $mainDir/BSD/Data
echo "Done.\nStart Annotation"
$mainDir/BSD/Script/Mut2Haplo.out $mainDir/BSD/Data/sample_mut_loc_time.tsv

echo "Done."
echo `date "+%Y-%m-%d %H:%M:%S"`