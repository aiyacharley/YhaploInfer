from __future__ import division
import os,sys
import csv
import gzip
from glob import glob
import argparse

def openFile(filename):
	if filename.endswith("gz"):
		handle = gzip.open(filename,"r")
	else:
		handle = open(filename,"r")
	return handle

def supportRatio(L):
	Ratio = 0
	for i in L[2].split(","):
		n = i.split("/")
		Ratio += float(int(n[0])/int(n[1]))
	return Ratio/int(L[0])
		
def supportRatioWeight(G,dictHaplo,dictT,mainSet):
#	G(O,O2,O2a,O2a1,O2a1b)
	a,b,c,d,e,f = [],[],[],[],[],[]
	scoreL = []
	for g in G.split(','):
		for i in dictT[g][2]:
			if i.split("/")[0] in mainSet:
				a.append(i)
			else:
				b.append(i)
		for j in dictT[g][0]:
			if j.split("/")[0] in mainSet:
				c.append(j)
			else:
				d.append(j)
		for k in dictHaplo[g]:
			if k.split("/")[0] in mainSet:
				e.append(k)
			else:
				f.append(k)
		if (len(e)+len(c))==0:
			aa = 1
		else:
			aa = len(e)+len(c)
		if (len(f)+len(d))==0:
	 		bb = 1	  
		else:
			bb = len(f)+len(d)
		score = 90*len(e)/aa+10*len(f)/bb
		scoreL.append(score)
	return sum(scoreL)/len(scoreL)

def main():		
	reader = csv.reader(mainSNP,delimiter='\t')
	mainSet = set()
	for rec in reader:
		for i in rec[1].split(", "):
			for j in i.split("/"):
				mainSet.add(j)
	
	set1 = set()
	for rec in csv.reader(ihandle1,delimiter='\t'):
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  16A0004401	
		if rec[0].startswith("#"):
			continue
		if rec[0]=="chrY" or rec[0]=="Y":
			ids = "|".join([rec[1],rec[3],rec[4]])
			set1.add(ids) 
	dictHaplo = {}
	for rec in csv.reader(ihandle2,delimiter='\t'):
	#2655180 G       A       O1b2    Page63/M176/SRY465
		ids = "|".join([rec[0],rec[1],rec[2]])
		if ids in set1:
			dictHaplo.setdefault(rec[3],set()).add(rec[4])
	#		print "\t".join(rec)
	dictT,dictG,dictRoot = {},{},{}
	setG = set()
	out2 = csv.writer(outF2,delimiter="\t")
	for rec in csv.reader(ihandle3,delimiter='\t'):
		dictT[rec[0]] = rec[1:]
	for i in sorted(dictHaplo.keys()):
		dictG[i] = [str(len(dictHaplo[i])),str(dictT[i][0]),str(dictT[i][2])]
		setG.add(i)
	#	group sample.N/1K.N/All.N sample.SNPs 1K.SNPs All.SNPs
	#	print i,str(len(dictHaplo[i]))+"/"+str(dictT[i][0])+"/"+str(dictT[i][2]),",".join(sorted(list(dictHaplo[i]))),dictT[i][1],dictT[i][3]
		out2.writerow([i,str(len(dictHaplo[i]))+"/"+str(dictT[i][0])+"/"+str(dictT[i][2]),",".join(sorted(list(dictHaplo[i]))),dictT[i][1],dictT[i][3]])
	dictR1,dictR2 = {},{}
	for rec in csv.reader(ihandle4,delimiter='\t'):
		if rec[0].startswith("Y"):
			dictRoot[rec[0].split(">")[-1]] = rec[0]
		Names = set(rec[0].split(">"))
		num = len(Names & setG)
		dictR1[rec[0]] = [num]
		dictR2[rec[0]] = [",".join(sorted(list(Names & setG)))]
	#	print num,",".join(sorted(list(Names & setG))),rec[0]
	out = csv.writer(outF,delimiter="\t")
	ind=0
	out3 = csv.writer(outP,delimiter="\t")
	dictPath,dictIMP = {},{}
	PathList = []
	for key,val in sorted(dictR1.items(), key=lambda d:d[1], reverse = True):
		if val[0]!=0:
			Slist = []
			for v in dictR2[key][0].split(","):
				Slist.append("/".join(dictG[v]))
	#		print val[0],",".join(sorted(dictR2[key])),",".join(Slist),key
			Ratio = supportRatio([val[0],",".join(dictR2[key]),",".join(Slist),key])
			score = supportRatioWeight(",".join(dictR2[key]),dictHaplo,dictT,mainSet)
			aa = dictR2[key][0].split(",")[-1]
			if ",".join(dictR2[key]) not in PathList:
				PathList.append(",".join(dictR2[key]))
				dictIMP[",".join(dictR2[key])] = float(score)*float(Ratio)*float(val[0])
			dictPath[",".join(dictR2[key])] = [score,Ratio,val[0],",".join(Slist),key[:key.index(aa)+len(aa)],",".join(dictR2[key])]
			out.writerow([score,Ratio,val[0],",".join(dictR2[key]),",".join(Slist),key])
			ind += 1
	for rec,val in sorted(dictIMP.items(), key=lambda d:d[1], reverse = True):
		if rec.split(",")[0]=="Y":
			continue
		if int(dictPath[rec][2])>3:
			out3.writerow([dictPath[rec][0],dictPath[rec][1],dictPath[rec][2],dictPath[rec][5],dictPath[rec][3],dictRoot[rec.split(",")[0]]+">"+dictPath[rec][4]])

if __name__=='__main__':
	parser = argparse.ArgumentParser(prog='python *.py',usage='%(prog)s -i input.vcf.gz -d /Path/to/output',description = 'infer Y Haplotype for 1K data',epilog = 'Created by WangCR. June 20, 2019')
	parser.add_argument('-i','--input',help='Input vcf format file.')
	parser.add_argument('-d','--outdir',default=".",help='Output file directory, defaulted current directory.')
	args = parser.parse_args()
	inputF = args.input
	outdir = args.outdir
	os.system("mkdir -p %s"%outdir)
	filename = inputF.split("/")[-1]
	ihandle1 = openFile(inputF)
	ihandle2 = openFile("Target1K_sites.txt.gz")
	ihandle3 = openFile("Target1K_HaploGroup.txt.gz")
	ihandle4 = openFile("TargetPath_ISOGG_2019_All.txt.gz")
	mainSNP = openFile("Target_mainSNPs.txt.gz")
	outF = gzip.open("%s/Result_%s"%(outdir,filename),"wb")
	outF2 = gzip.open("%s/Result_sites_%s"%(outdir,filename),"wb")
	outP = gzip.open("%s/Result_Path_%s"%(outdir,filename),"wb")
	main()
	ihandle1.close()
	ihandle2.close()
	ihandle3.close()
	ihandle4.close()
	outF.close()
	outF2.close()
	outP.close()
	mainSNP.close()
