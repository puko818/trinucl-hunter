#!/bin/bash

#download the repeats data

#http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
#http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz


UPSTREAM=500
DOWNSTREAM=500

trap 'echo "error"; do_cleanup failed; exit' ERR
trap 'echo "received signal to stop"; do_cleanup interrupted; exit' SIGQUIT SIGTERM SIGINT

do_cleanup () { echo "$1 $(date)" >> script_log; }

mkdir -p "RAW_DATA"

download_repeats()
{
echo "Downloading simple repeats..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" | gunzip -c > "RAW_DATA/"simpleRepeat.txt

echo "Downloading repeatmasker repeats..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" | gunzip -c > "RAW_DATA/"rmRepeat.txt

echo "Downloading microsattelite repeats..."
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/microsat.txt.gz" | gunzip -c > "RAW_DATA/"microsatRepeat.txt
}



#CREATE TABLE `simpleRepeat` (
#`bin` smallint(5) unsigned NOT NULL,		$1
#`chrom` varchar(255) NOT NULL,       		$2
#`chromStart` int(10) unsigned NOT NULL,	$3
#`chromEnd` int(10) unsigned NOT NULL,		$4
#`name` varchar(255) NOT NULL,			$5
#`period` int(10) unsigned NOT NULL,		$6
#`copyNum` float NOT NULL,			$7
#`consensusSize` int(10) unsigned NOT NULL,	$8
#`perMatch` int(10) unsigned NOT NULL,		$9
#`perIndel` int(10) unsigned NOT NULL,		$10
#`score` int(10) unsigned NOT NULL,		$11
#`A` int(10) unsigned NOT NULL,			$12
#`C` int(10) unsigned NOT NULL,			$13
#`G` int(10) unsigned NOT NULL,			$14
#`T` int(10) unsigned NOT NULL,			$15
#`entropy` float NOT NULL,			$16
#`sequence` longblob NOT NULL,			$17
#KEY `chrom` (`chrom`(16),`bin`),
#KEY `chrom_2` (`chrom`(16),`chromStart`)

#cat "RAW_DATA"/simpleRepeat.txt | egrep -v "A|T" > "RAW_DATA"/simpleRepeatGC.txt
#cat "RAW_DATA"/simpleRepeat.txt | egrep -v "G|C" > "RAW_DATA"/simpleRepeatAT.txt



#CREATE TABLE `rmsk` (
#  `bin` smallint(5) unsigned NOT NULL,		$1
#  `swScore` int(10) unsigned NOT NULL,		$2
#  `milliDiv` int(10) unsigned NOT NULL,	$3
#  `milliDel` int(10) unsigned NOT NULL,	$4
#  `milliIns` int(10) unsigned NOT NULL,	$5
#  `genoName` varchar(255) NOT NULL,		$6
#  `genoStart` int(10) unsigned NOT NULL,	$7
#  `genoEnd` int(10) unsigned NOT NULL,		$8
#  `genoLeft` int(11) NOT NULL,			$9
#  `strand` char(1) NOT NULL,			$10
#  `repName` varchar(255) NOT NULL,		$11
#  `repClass` varchar(255) NOT NULL,		$12
#  `repFamily` varchar(255) NOT NULL,		$13
#  `repStart` int(11) NOT NULL,			$14
#  `repEnd` int(11) NOT NULL,			$15
#  `repLeft` int(11) NOT NULL,			$16
#  `id` char(1) NOT NULL,			$17
#  KEY `genoName` (`genoName`(14),`bin`)


#CREATE TABLE `microsat` (
#  `bin` smallint(5) unsigned NOT NULL,
#  `chrom` varchar(255) NOT NULL,
#  `chromStart` int(10) unsigned NOT NULL,
#  `chromEnd` int(10) unsigned NOT NULL,
#  `name` varchar(255) NOT NULL,
#  KEY `name` (`name`(16)),
#  KEY `chrom` (`chrom`(14),`bin`)
#) ENGINE=MyISAM DEFAULT CHARSET=latin1;
#/*!40101 SET character_set_client = @saved_cs_client */;


#cat "RAW_DATA"/rmRepeat.txt
#cat "RAW_DATA"/rmRepeat.txt


ATNAME="AT"
GCNAME="GC"

DBFOLDER="DB"
FINAL="FINAL"

mkdir -p $FINAL
mkdir -p $DBFOLDER

create_database()
{
echo "Creating database of GCrich and ATrich repeats .."
{ cat "RAW_DATA"/rmRepeat.txt | grep Simple_repeat | ./processRmsk.awk; cat "RAW_DATA"/simpleRepeat.txt | ./processSimple.awk; cat "RAW_DATA/"microsatRepeat.txt | ./processMicrosat.awk; } |sort   -t$'\t' -k1 -k2,3n -k7,7n  | sort -t$'\t' -u  -k1,1 -k2,3n  |tee >(egrep -v "A|T"  | ./remove_overlapping_intervals.py | bedtools sort -i - > $DBFOLDER/$GCNAME.ALL.ALL.txt) |  egrep -v "G|C" | ./remove_overlapping_intervals.py | bedtools sort -i - >  $DBFOLDER/$ATNAME.ALL.ALL.txt

echo "Removing overlapping repeats..."
for i in $ATNAME $GCNAME
do
 cat $DBFOLDER/$i".ALL.ALL.txt" | tee >(awk -F$"\t" '$4==3 {print}' > $DBFOLDER/$i".3.ALL.txt") |\
                      awk -F$"\t" '$4==3 && $5>50 {print}' > $DBFOLDER/$i".3.50.txt"
done
}

download_annotations_hg38() {
echo "Downloading data for the overlapping analysis..."

echo "Downloading introns annotation from hg38..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz" |gunzip -c  |\
   awk 'BEGIN{OFS="\t"}{ exonCount=int($8);split($9,exonStarts,"[,]"); split($10,exonEnds,"[,]"); for(i=1;i<exonCount;i++) {printf("%s\t%d\t%d\t%s\n",$2,exonEnds[i],exonStarts[i+1],$1);}}' |\
   bedtools sort -i - \
   > "RAW_DATA"/introns.bed

echo "Downloading exons annotation from hg38..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz" |gunzip -c  |\
   awk 'BEGIN{OFS="\t"}{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s\t%d\t%d\t%s\n",$2,S[i],E[i],$1);} }' |\
   bedtools sort -i - \
   > "RAW_DATA"/exons.bed


echo "Downloading promoters annotation from hg38..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz" |gunzip -c |\
   awk 'BEGIN {OFS="\t"}{ if ($3 == "+") { print $2,($4-'$UPSTREAM' < 0 ? 0 :$4-'$UPSTREAM'), $4+'$DOWNSTREAM',$1 } else if ($6 == "-") { print $2,($5-'$UPSTREAM' < 0 ? 0 :$5-'$UPSTREAM'), $5+'$DOWNSTREAM',$1 }}' |\
   bedtools sort -i - \
   > "RAW_DATA"/promoters.bed


echo "Downloading cytogenetic bands annotation from hg38..."
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz" | gunzip -c |\
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' |\
    bedtools sort -i - \
    > "RAW_DATA"/cyto.bed

echo "Downloading transposable elements annotation from hg38"
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" | gunzip -c  |\
   grep -v "Simple_repeat" | \
   awk 'BEGIN {OFS="\t"} {print $6,$7,$8,$11"_"$12}'  |\
   bedtools sort -i - \
   > "RAW_DATA"/te.bed

echo "Downloading CpG islands annotation from hg38"
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz" | gunzip -c | awk -F"\t" 'BEGIN {OFS="\t"}{gsub(" ", "", $5);print $2,$3,$4,$5}' | bedtools sort -i - > "RAW_DATA"/cpg.bed

echo "Downloading SNP annotation from hg38"
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz" | gunzip -c | awk -F"\t" 'BEGIN {OFS="\t"}{id=$5"_"$12; if ($3==$4) print $2,$3,$4+1,id; else print $2,$3,$4,id; }' | bedtools sort -i - > "RAW_DATA"/snp.bed

echo "Download DGV struct from hg38"
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/dgvMerged.txt.gz" | gunzip -c | awk  -F"\t" 'BEGIN {OFS="\t"}{id=$5"_"$11; if ($3==$4) print $2,$3,$4+1,id; else print $2,$3,$4,id;}' | bedtools sort -i - > "RAW_DATA"/dgv.bed

}


download_annotations_hg19()
{
curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRegTfbsClusteredV3.txt.gz" | gunzip -c | awk -F"\t" 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}' | liftOver stdin  hg19ToHg38.over.chain stdout wgEncodeRegTfbsClusteredV3.unmapped > "RAW_DATA"/wgEncodeRegTfbsClusteredV3.bed

#curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRegTxnCaltechRnaSeqHelas3R2x75Il200SigPooled.txt.gz" | gunzip -c

curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cosmic.txt.gz" | gunzip -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5}' |  liftOver stdin  hg19ToHg38.over.chain stdout cosmic.unmapped > "RAW_DATA"/cosmic.bed

}


create_header()
{
echo -ne "chrom\t"
echo -ne "chromStart\t"
echo -ne "chromEnd\t"
echo -ne "period\t"
echo -ne "copyNum\t"
echo -ne "totalLen\t"
echo -ne "perMatch\t"
echo -ne "sequence\t"
echo -ne "source\t"
echo -ne "cytoBand\t"
echo -ne "cytoBand:L\t"
echo -ne "CpG\t"
echo -ne "CpG:L\t"
echo -ne "promoter\t"
echo -ne "promoter:L\t"
echo -ne "exon\t"
echo -ne "exon:L\t"
echo -ne "intron\t"
echo -ne "intron:L\t"
echo -ne "TE\t"
echo -ne "TE:L\t"
echo -ne "TFBS\t"
echo -ne "TFBS:L\t"
echo -ne "cosmic\t"
echo -ne "cosmic:L\t"
echo -ne "dbSNP\t"
echo -ne "dbSNP:L\t"
echo -ne "DGV\t"
echo "DGV:L"
}

decorate_with_annotations()
{
bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size $1 "RAW_DATA"/cyto.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/cpg.bed  |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/promoters.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/exons.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/introns.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/te.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/wgEncodeRegTfbsClusteredV3.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/cosmic.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/snp.bed |\
                     bedmap  --bp-ovr 1  --echo --echo-map-id --echo-overlap-size - "RAW_DATA"/dgv.bed
}

#create_trinucleotide_datasets()
#{
#  ar=("$@")
#  for ((i=0;i<${#ar[@]};++i)); 
#  do
#    #echo ${ar[i]}
#    cat ${ar[i]} | awk -F$'\t' '$4==3 && $5>50 {print}'   > ${ar[i]%*.txt}".trinucl50.txt"
#  done
#}


analyse_overlaps()
{

  repeat=`basename $1 | cut -d. -f1`
  period=`basename $1 | cut -d. -f2`
  copyn=`basename  $1 | cut -d. -f3`
  genomesize=`curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$2 }END{print SUM}'`
  repeatspan=`sortBed -i $1| mergeBed -i - | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'` 
  featurespan=`sortBed -i $2 | mergeBed -i - | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'`
  overlap=`bedmap  --bp-ovr 1   --echo-overlap-size $1 <(sortBed -i $2 | mergeBed -i - | sortBed -i -) |  awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$1 }END{print SUM}'`
  echo -e "$repeat\t$period\t$copyn\t$3\t$repeatspan\t$featurespan\t$overlap\t$genomesize"
}

analyse_overlaps_header()
{
 echo -e "repeat\tperiod\tcopyn\tfeature\trepeatspan\tfeaturespan\toverlap\tgenomesize"
}

analyse_overlaps_for_plots()
{
analyse_overlaps_header

features=( "exons.bed" "introns.bed" "te.bed" "promoters.bed" "cosmic.bed" "wgEncodeRegTfbsClusteredV3.bed" "snp.bed" "dgv.bed"  )
feature_names=("Exons" "Introns" "TEs" "Promoters" "COSMIC" "TFBS" "SNPs" "Structural variations")
###dbs=($GCNAME $ATNAME)
###copyn
###period

###dbnames=("GC-rich" "AT-rich")
###dbs_filtered=("GCrich.no_overlaps.trinucl50.txt" "ATrich.no_overlaps.trinucl50.txt")
###dbnames_filtered=("GC-rich-3nt" "AT-rich-3nt")

###create_trinucleotide_datasets "${dbs[@]}"


for j in $DBFOLDER/*;
 do
  for ((i=0;i<${#features[@]};++i));
   do
     analyse_overlaps $j  RAW_DATA/${features[i]} ${feature_names[i]}
   done
 done

}


#####################################
#THIS IS THE ENTRY TO THE PROGRAMME


download_repeats

create_database



download_annotations_hg38

download_annotations_hg19



for i in $GCNAME $ATNAME
do
  echo "Decorating "$i"..."
  create_header > $FINAL/$i"_with_features.txt"
  #decorate_with_annotations $DBFOLDER/$i".txt" | sed 's/|/\t/g' |  awk -F$'\t' '$4==3' | sort   -t$'\t'  -k6,6nr -k5,5nr >> $DBFOLDER/$i"_with_features.txt"
  decorate_with_annotations $DBFOLDER/$i".ALL.ALL.txt" | sed 's/|/\t/g' |  sort   -t$'\t'  -k6,6nr -k5,5nr >> $FINAL/$i"_with_features.txt"
done

#echo "Generating data for plots ..."

#analyse_overlaps_for_plots > $FINAL/"EXP_OBS.txt"

#to get only period 3 sorted by total len and copynum:
#cat GCrich.txt | ./remove_overlapping_intervals.py |  awk -F$'\t' '$4==3' | sort   -t$'\t'  -k6,6nr -k5,5nr
#cat ATrich.txt | ./remove_overlapping_intervals.py |  awk -F$'\t' '$4==3' | sort   -t$'\t'  -k6,6nr -k5,5nr

#cat "RAW_DATA"/simpleRepeat.txt | ./processSimple.awk
#cat "RAW_DATA/"microsatRepeat.txt | ./processMicrosat.awk





