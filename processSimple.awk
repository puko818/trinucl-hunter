#!/bin/awk -f

#CREATE TABLE `simpleRepeat` (
#`bin` smallint(5) unsigned NOT NULL,           $1
#`chrom` varchar(255) NOT NULL,                 $2
#`chromStart` int(10) unsigned NOT NULL,        $3
#`chromEnd` int(10) unsigned NOT NULL,          $4
#`name` varchar(255) NOT NULL,                  $5
#`period` int(10) unsigned NOT NULL,            $6
#`copyNum` float NOT NULL,                      $7
#`consensusSize` int(10) unsigned NOT NULL,     $8
#`perMatch` int(10) unsigned NOT NULL,          $9
#`perIndel` int(10) unsigned NOT NULL,          $10
#`score` int(10) unsigned NOT NULL,             $11
#`A` int(10) unsigned NOT NULL,                 $12
#`C` int(10) unsigned NOT NULL,                 $13
#`G` int(10) unsigned NOT NULL,                 $14
#`T` int(10) unsigned NOT NULL,                 $15
#`entropy` float NOT NULL,                      $16
#`sequence` longblob NOT NULL,                  $17
#KEY `chrom` (`chrom`(16),`bin`),
#KEY `chrom_2` (`chrom`(16),`chromStart`)




BEGIN{OFS="\t"} 
{
  print $2,$3,$4,$6,$7,$6*$7,$9,$17,"trf"
}



