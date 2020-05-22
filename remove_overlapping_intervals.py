#!/usr/bin/env python
#DESCRIPTION
# script takes chrom|chromstart|chromEnd|period|copyNum|totalLen|perMatch|seq|program format
# and in case of overlapping records reports the one with higher perMatch
#
#


import fileinput

stack=[]

def isOverlap(interval1,interval2):
  retval=0
  if interval1[0]==interval2[0]:
    latest_start=max(interval1[1],interval2[1])
    earliest_end=min(interval1[2],interval2[2])
    delta=earliest_end-latest_start
    if delta>0:
      retval=delta+1
  return retval
      




for line in fileinput.input():
    actl=line.strip().split("\t")
    if not stack:
      stack.append(actl)
    else:
      interval1=(str(stack[-1][0]),int(stack[-1][1]),int(stack[-1][2]))
      interval2=(str(actl[0]),int(actl[1]),int(actl[2]))
      if isOverlap(interval1,interval2)>0:
        score1=float(stack[-1][5])
        score2=float(actl[5])
        period1=int(stack[-1][3])
        period2=int(actl[3])
        if period1==period2:#if the same period decide entirely based on the score
          if score2>=score1:
             stack[-1]=actl
        elif period2==3:#if the period of the new record is 3 and the previos record is not 3, then 3 has priority
          stack[-1]=actl
        elif period1!=3:#if they are not equal but neither of them is 3 decide based on the score
          if score2>=score1:
            stack[-1]=actl
      else:
        print "\t".join(stack[-1])
        del(stack[-1])
        stack.append(actl)

if stack:
 print "\t".join(stack[-1])
        
    
     

   
