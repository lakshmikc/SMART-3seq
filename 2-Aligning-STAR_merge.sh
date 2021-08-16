
#!/bin/sh


for name in PR1643*; do

echo $name

python umi_homopolymer.py ./$name/L001/*_R1_*.fastq ./$name/L001/$name\_L001.fastq
python umi_homopolymer.py ./$name/L002/*_R1_*.fastq ./$name/L002/$name\_L002.fastq
python umi_homopolymer.py ./$name/L003/*_R1_*.fastq ./$name/L003/$name\_L003.fastq
python umi_homopolymer.py ./$name/L004/*_R1_*.fastq ./$name/L004/$name\_L004.fastq

STAR --runThreadN 60 --quantMode GeneCounts --genomeDir PIG_Data/Pig_genome/star_index --readFilesIn  ./$name/L001/${name}_L001.fastq   --outFileNamePrefix ./$name/L001/${name}_L001_star  --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --outFilterMatchNmin 0 --outFilterMismatchNmax 999 -clip3pAdapterMMP 0.2 -clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --twopassMode Basic
STAR --runThreadN 60 --quantMode GeneCounts --genomeDir PIG_Data/Pig_genome/star_index --readFilesIn  ./$name/L002/${name}_L002.fastq   --outFileNamePrefix ./$name/L002/${name}_L002_star  --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --outFilterMatchNmin 0 --outFilterMismatchNmax 999 -clip3pAdapterMMP 0.2 -clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --twopassMode Basic
STAR --runThreadN 60 --quantMode GeneCounts --genomeDir PIG_Data/Pig_genome/star_index --readFilesIn  ./$name/L003/${name}_L003.fastq   --outFileNamePrefix ./$name/L003/${name}_L003_star  --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --outFilterMatchNmin 0 --outFilterMismatchNmax 999 -clip3pAdapterMMP 0.2 -clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --twopassMode Basic
STAR --runThreadN 60 --quantMode GeneCounts --genomeDir Pig_genome/star_index --readFilesIn  ./$name/L004/${name}_L004.fastq   --outFileNamePrefix ./$name/L004/${name}_L004_star  --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --outFilterMatchNmin 0 --outFilterMismatchNmax 999 -clip3pAdapterMMP 0.2 -clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --twopassMode Basic


java -jar ~/Downloads/picard-2.18.14.jar MergeSamFiles I=./$name/L001/$name\_L001_star_mmAligned.out.sam I=./$name/L002/$name\_L002_star_mmAligned.out.sam I=./$name/L003/$name\_L003_star_mmAligned.out.sam I=./$name/L004/$name\_L004_star_mmAligned.out.sam O=/home/lck003/SPARC/PIG_Data/PR1643-97/STAR/$name\.sam

done
