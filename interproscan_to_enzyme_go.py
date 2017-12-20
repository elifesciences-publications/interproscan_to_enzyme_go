#! /usr/bin/python
import time
import argparse
import re
import goatools #https://github.com/tanghaibao/goatools,
#!wget http://purl.obolibrary.org/obo/go/go-basic.obo",


parser = argparse.ArgumentParser(description='Parse the InterProScan tsv file and output all those with GO term that are children of the given GO parent term')
parser.add_argument('interproscan_tsv',help="file with the uniprot IDs to lookup gene ontology for")
parser.add_argument('-o',help="Output file")
parser.add_argument('--go_term',default="GO:0003824",help="GO term to look up parents for. Default=GO:0003824 - enzymatic activity")
parser.add_argument('--override_file',default=None,help="CSV file with column1:string to search for and column2:GO_term to assign")

args = parser.parse_args()

print args

print "Loading gene ontology directed acyclic graph..."
go_dag = goatools.obo_parser.GODag("go-basic.obo")
print "Finished loading."

override = dict()
if args.override_file != None:
	handle = open(args.override_file)
	for line in handle:
		splitline = line.split(",")
		string = splitline[0].strip()
		GO = splitline[1].strip()
		if string not in override.keys():
			override[string] = [GO]
		else:
			override[string].append(GO)
	print "Loaded",len(override),"annotation overrides."


def GO_has_parent(child_GO_term,parent_GO_term):
	child_GO_obj = go_dag[child_GO_term]
	if parent_GO_term in child_GO_obj.get_all_parents():
		return True

	elif child_GO_term == parent_GO_term:
		return True
	else:
		return False 

if args.override_file != None:
	override_regex_string = "|".join(override.keys())
	override_regex = re.compile("|".join(override.keys()))
go_regex = re.compile("GO:[0-9]{7}|GO:[0-9]{8}")

gene_to_GO_dict = dict()

j=0
handle = open(args.interproscan_tsv,"rU")
for line in handle.readlines():
	if line[0] == "#":
		continue
	splitline = line.split("\t")
	gene_id = splitline[0].strip()
	matches = go_regex.findall(line)
	if len(matches) != 0:
		for m in matches:
			if gene_id not in gene_to_GO_dict.keys():
				gene_to_GO_dict[gene_id] = set()
				gene_to_GO_dict[gene_id].add(m)
			else:
				gene_to_GO_dict[gene_id].add(m)
	if args.override_file != None:
		override_matches = override_regex.findall(line)
	else:
		override_matches = []
	if len(override_matches) != 0:
		for m in override_matches:
			j+=1
                        if gene_id not in gene_to_GO_dict.keys():
                                gene_to_GO_dict[gene_id] = set()
				for g in override[m]:
                                	gene_to_GO_dict[gene_id].add(g)
                        else:
				for g in override[m]:
                                	gene_to_GO_dict[gene_id].add(g)
out_handle = open(args.o,"wb")
i=0
completed=set()
print "Looking up genes whose GO terms are children of:",args.go_term
for id in gene_to_GO_dict.keys():
	for go in gene_to_GO_dict[id]:
		if GO_has_parent(go,args.go_term) and id not in completed:
			out_handle.write(id+"\t"+args.go_term+"\t"+",".join(gene_to_GO_dict[id])+"\n")
			i+=1
			completed.add(id)
			continue
out_handle.close()
if args.override_file != None:
	print "Made",j,"override matches."
print "Wrote",i,"gene ids."

