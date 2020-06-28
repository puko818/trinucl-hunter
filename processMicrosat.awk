#!/bin/awk -f

#CREATE TABLE `microsat` (
#  `bin` smallint(5) unsigned NOT NULL, 	$1
#  `chrom` varchar(255) NOT NULL,  		$2
#  `chromStart` int(10) unsigned NOT NULL,	$3
#  `chromEnd` int(10) unsigned NOT NULL,	$4
#  `name` varchar(255) NOT NULL,		$5  [0-9]+x[ATGC]+
#  KEY `name` (`name`(16)),			$6  nnx[ATGC]
#  KEY `chrom` (`chrom`(14),`bin`)		$7
#) ENGINE=MyISAM DEFAULT CHARSET=latin1;


BEGIN{OFS="\t"} 
{
  match($5, /([0-9]+)x([ATGC]+)/,a) 
  period=length(a[2]); 
  copynum=a[1]
  percmatch=100; 
  totallen=period*copynum
  print $2,$3,$4,period,copynum,totallen,percmatch,a[2],"microsat"
}



