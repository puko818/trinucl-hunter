#!/bin/awk -f
BEGIN {}  # Begin section
{}        # Loop section
END{}     # End section

BEGIN{OFS="\t"} 
{
  match($11, /([ATGC]+)/,a)
  #totallen=$15-$14; 
  totallen=$8-$7
  period=length(a[1]); 
  copynum=totallen/period; 
  percmatch=(1-(($3+$4+$5)/1000))*100; 
  print $6,$7,$8,period,copynum,totallen,percmatch,a[1],"rpmsk"
}



