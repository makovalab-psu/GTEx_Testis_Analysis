#!/usr/bin/python
import sys


Test_input="/galaxy/home/rxv923/src/yhaplo/data/1000Y.subset.genos.txt"
Inp=open(Test_input, "rU")

locations=Inp.readline()
col=locations.split("\t")
del(col[0])
results = map(int, col)
interest=dict()
for i in results:
	if i in interest:
		interest[i]+=1
	else:
		interest[i]=0


value=dict()
Input=sys.argv[1]###"SRR2157516_Ychr_hg19_chrY_2655100_28771000.pileup"

#Parse the pileup file. Since Y is haploid, for each position we obtain the frequency of each nucleotide based on the read depth. If the frequency is above 0.7 then we assign that nucleotide and if none of them have high frequency we return the reference at that position. This is a naive implementation and can be improved in the future.
with open(Input, "rU") as f:
	for line in f:
		if line[0]=="#":
			continue
		else:
			A,T,G,C=0,0,0,0
			col=line.split("\t")
			chr=col[0]
			position=col[1]
			ref=col[2]
			depth=col[3]
			pile=col[4].upper()
			if float(depth) > 0 :
				A=pile.count('A')
				T=pile.count('T')
				G=pile.count('G')
				C=pile.count('C')
				Ins=pile.count('+')
				Del=pile.count('-')
				if (Ins+Del)/float(depth) <= 0.2:
					if A/float(depth) > 0.7:
						value[position]="A"
					elif T/float(depth) > 0.7:
						value[position]="T"
					elif G/float(depth) > 0.7:
						value[position]="G"
					elif C/float(depth) > 0.7:
						value[position]="C"
					else:
						value[position]=ref.upper()
				else:
					value[position]="-"
			else:
				value[position]="-"

########Create thegenos file by reading the test file and parsing the nucleotide in the locations defined in the test file.
header=list("ID")
row=list("Sample")
for i in interest.keys():
	if str(i) in value:
		header.append("\t"+str(i))
		row.append("\t"+value[str(i)])
	else:
		print i


Out=sys.argv[2]##'SRR2157516_Ychr_hg19.genos.txt'
f = open(Out, 'w')
f.write("".join(header)+'\n')  
f.write("".join(row)+'\n')  
f.close()#!/usr/bin/python

