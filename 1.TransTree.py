import os,sys
import csv,re
import itertools
import networkx as nx
from glob import glob

def createEdge(edgeList,net):
	L = len(edgeList)
	edge = edgeList[::-1]
	for indx,i in enumerate(edge):
		for indj,j in enumerate(edge[indx+1:]):
			if int(i)>int(j):
				net.add_edge(L-indj-indx-2,L-indx-1)
				break
	return net

def processISOGG_TreeFile(f): 
	net = nx.DiGraph()
	edgeList = []
	ind = 0
	for line in f:
		line = line.rstrip()
		if len(line)==0:
			continue
		if "~" in line:
			continue
		if ind==0:
			line = str(ind)+"\t0\t"+line.strip()
		else:
			l = [(k, len(list(g))) for k, g in itertools.groupby(line)]
			line = str(ind)+"\t"+str(l[0][1])+"\t"+line.strip()
		line = line.split("\t")
		if len(line)==3:
			line = line.append("NA")
		if line==None:
			continue
		if "or" in line[2]:
			Nname = line[2].split(" or")[0].strip().split(" ")[0]
		else:
			Nname = line[2].strip().split(" ")[0]
		net.add_node(ind, name=Nname, DeriveSNP=line[3].replace(", ",","))
		edgeList.append(line[1])
		ind += 1
	return edgeList,net

def main():
	edgeList,net = processISOGG_TreeFile(f)
#	print edgeList
	net = createEdge(edgeList,net)
#	print net.adj
#	print net._node
	Node1,Node2 = set(),set()
	for n,nbrs in net.adjacency():
		for nbr,attr in nbrs.items():
#			out.writerow([n,nbr,net.node[n]["name"],net.node[nbr]["name"],net.node[nbr]["DeriveSNP"]])
			Node1.add(n)
			Node1.add(nbr)
			Node2.add(n)
	leaveNode = Node1 - Node2
	for i in sorted(leaveNode):
		path = nx.dijkstra_path(net,source=0,target=i)
		plist = []
		for p in path:
			plist.append(net.node[p]["name"])
#		print "->".join(plist)
		out2.writerow([">".join(plist)])

if __name__ == '__main__':
	files = glob("ISOGG2019/ISOGG_2019_*")
	for inFile in files:
#		inFile = sys.argv[1]
		f = open(inFile, "r")
		fname = inFile.split("/")[-1]
#		ohandle = open("Tree_%s"%sys.argv[1],"wb")
#		out = csv.writer(ohandle,delimiter="\t")
		ohandle2 = open("Path_%s"%fname,"wb")
		out2 = csv.writer(ohandle2,delimiter="\t")
		main()
		f.close()
#		ohandle.close()
		ohandle2.close()
	os.system("cat Path* > TargetPath_ISOGG_2019_All.txt")
	os.system("rm Path_*")

