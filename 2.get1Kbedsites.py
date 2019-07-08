import os,sys
import csv
from glob import glob

handle1 = open("1Kgene_capture_targets.chrY.bed","rU") # 1Kgene_capture_targets.chrY.bed
Ylist = [0]*28818010  #28818008
for rec in csv.reader(handle1,delimiter="\t"):
	for i in range(int(rec[1]),int(rec[2])+1):
		Ylist[i] = 1
infiles = glob("mutation_ISOGG2019_*.txt")
dictTargetSNP1 = {}
dictTargetSNP2 = {}
dictPos = {}
dictPosAll = {}
dictAll,dictAll2 = {},{}
for inF in infiles:
	handle2 = open(inF,"rU") # mutation_ISOGG2019_HaploA.txt
	header = handle2.next().strip().split("\t")
	for rec in csv.reader(handle2,delimiter="\t"):
#	Name;Haplogroup;Other Names;rs #;Build 37 #;Build 38 #;Mutation info
		if "~" in rec[1]:
			continue
		if len(rec[6])==0:
			continue
		if rec[4]=='None':
			continue
		if "ins" in rec[6]:
			continue
		position = rec[4].split("..")[0]
		if rec[6].endswith("->del"):
			bases = rec[6].replace("del","-")
		elif rec[6].startswith("del->"):
			bases = rec[6].replace("del","-")
		else:	
			bases = rec[6]
		siteID = "|".join([position]+bases.split('->'))
		dictAll.setdefault(siteID,set()).add(rec[0])
		dictAll2.setdefault(siteID,set()).add(rec[1].split(" ")[0])
		dictPosAll[int(position)] = siteID
		if Ylist[int(position)] == 1:
			dictPos[int(position)] = siteID
			dictTargetSNP1.setdefault(siteID,set()).add(rec[0])
			dictTargetSNP2.setdefault(siteID,set()).add(rec[1].split(" ")[0])
#			dictTargetSNP2.setdefault(siteID,set()).add(rec[1])
	handle2.close()	
ohandle = open("Target1K_sites.txt","wb")
out = csv.writer(ohandle,delimiter='\t')
dictHaplo = {}
dictHaploAll = {}
for pos in sorted(dictPosAll.keys()):
	ids = dictPosAll[pos]
	V1 = list(dictAll2[ids])
	V2all = "/".join(sorted(list(dictAll[ids])))
	for v1 in V1:
		if "or" in v1:
			name = v1.split(" or")[0]
		else:
			name = v1
		dictHaploAll.setdefault(name,[]).append(V2all)

for pos in sorted(dictPos.keys()):
	ids = dictPos[pos]
	V1 = list(dictTargetSNP2[ids])
	V2 = "/".join(sorted(list(dictTargetSNP1[ids])))
	for v1 in V1:
		if "or" in v1:
			name = v1.split(" or")[0]
		else:
			name = v1
		out.writerow(ids.split("|")+[name,V2])
		dictHaplo.setdefault(name,[]).append(V2)
ohandle2 = open("Target1K_HaploGroup.txt","wb")
out2 = csv.writer(ohandle2,delimiter="\t")
for h in sorted(dictHaplo.keys()):
#	print h,len(dictHaplo[h]),len(dictHaploAll[h])
	out2.writerow([h,len(dictHaplo[h]),",".join(dictHaplo[h]),len(dictHaploAll[h]),",".join(dictHaploAll[h])])

handle1.close()
