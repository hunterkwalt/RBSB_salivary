#3/15/2023
#this is a script to explore OMA analysis with Pyham
import pyham
import re
import pandas as pd

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

import os
import sys

#take the newick and the orthoxml as arguments from the command line
nwk_path = sys.argv[1]
tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
orthoxml_path = sys.argv[2]

ham_analysis = pyham.Ham(tree_str, orthoxml_path)

treeprofile = ham_analysis.create_tree_profile(outfile="Output/tp.html")


# Get the genome of interest
pg = ham_analysis.get_extant_genome_by_name("Piezodorus_guildinii")
pg_nv_hh = ham_analysis.get_ancestral_genome_by_name("Piezodorus_guildinii/Nezara_viridula/Halyomorpha_halys")
pg_nv_hh_eu = ham_analysis.get_ancestral_genome_by_name("Piezodorus_guildinii/Nezara_viridula/Halyomorpha_halys/Euschistus_heros")
nv = ham_analysis.get_extant_genome_by_name("Nezara_viridula")
hh = ham_analysis.get_extant_genome_by_name("Halyomorpha_halys")
eh = ham_analysis.get_extant_genome_by_name("Euschistus_heros")

# Instanciate the gene mapping !
pentatomids = ham_analysis.compare_genomes_vertically(pg,pg_nv_hh) # The order doesn't matter!
pg_nv = ham_analysis.compare_genomes_lateral(pg,nv)
pg_hh = ham_analysis.compare_genomes_lateral(pg,hh)
pg_eh = ham_analysis.compare_genomes_lateral(pg,eh)

# The duplicated genes (that have duplicated) 
#print("HOG at vertebrates -> list of descendants gene in human")
#print(pentatomids.get_lost())
#print("\n")

# The gained genes (that emerged in between)
#print("List of human gene")
#print(pentatomids.get_gained())
#print("\n")

# The lost genes (that been lost in between) 
#print("HOG at vertebrates that are lost")
#print(pentatomids.get_lost())
#print("\n")


#get genes that have been gained since stink bug divergence
vertical_stink_bugs_pg = ham_analysis.compare_genomes_vertically(pg, pg_nv_hh)



# The duplicated genes (that have duplicated) 
pg_dup = vertical_stink_bugs_pg.get_duplicated()
dup_of = open("duplicated_pg_specific", "w")
#for g in pg_dup:
	#print(g)
#	gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(g))
#	print(gene_num)
#	#gene_name = ham_analysis.get_gene_by_id(gene_num[0])
#	dup_dict = gene_name.get_dict_xref()
#	dup_genes = dup_dict['protId']
#	dup_of.write(dup_genes+"\n")


for l in vertical_stink_bugs_pg.get_duplicated().values():
	for element in l:
		g = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))[0]
		dup_genes = ham_analysis.get_gene_by_id(g)
		dup_dict = dup_genes.get_dict_xref()
		dup_genes = dup_dict['protId']
		dup_of.write(dup_genes+"\n")
		
print("\n")


# The gained genes (that emerged in between)
pg_gained = vertical_stink_bugs_pg.get_gained()

#make a file of the genes names that were gained in the RBSB specific lineage
gained_of = open("gained_pg_specific", "w")
for g in pg_gained:
	#print(g)
	gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(g))
	gene_name = ham_analysis.get_gene_by_id(gene_num[0])
	gained_dict = gene_name.get_dict_xref()
	gained_genes = gained_dict['protId']
	gained_of.write(gained_genes+"\n")

#get genes that were lost	
#I will have to do a little more work on this to extract the names of the HOGs that the genes came from, but for now, I am going to grab them after I get the gene list.
pg_lost = vertical_stink_bugs_pg.get_lost()
lost_of = open("lost_pg_specific", "w")
for h in pg_lost:
	lost_genes = h.get_all_descendant_genes()
	for element in lost_genes:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		lost_dict = gene_name.get_dict_xref()
		lost_genes = lost_dict['protId']
		lost_of.write(lost_genes+"\n")
	


#nv_lost_pg = pg_nv.get_lost()
#hh_lost_pg = pg_hh.get_lost()

'''nv_lost = open("nv_vs_pg_lost", "w")
for h in nv_lost_pg:
	lost_genes = h.get_all_descendant_genes()
	for element in lost_genes:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		lost_dict = gene_name.get_dict_xref()
		lost_genes = lost_dict['protId']
		lost_of.write(lost_genes+"\n")'''
		
# The lost genes (that been lost in between) 
nv_hogs_lost = []
print("LOST GENES")
for hog, genomes in pg_nv.get_lost().items():
    print("\t- {} at pg/nv/hh has been lost in ".format(hog))
    if len(genomes) == 1:
    	for g in genomes:
    		if format("{}".format(g)) == "Nezara_viridula":
    			nv_hogs_lost.append(hog)
print("\n")

nv_lost = open("nv_vs_pg_lost", "w")
for h in nv_hogs_lost:
	lost_genes = h.get_all_descendant_genes()
	for element in lost_genes:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		lost_dict = gene_name.get_dict_xref()
		lost_genes = lost_dict['protId']
		nv_lost.write(lost_genes+"\n")


#now get the genes lost in halyomorpha halys
hh_hogs_lost = []
print("LOST GENES")
for hog, genomes in pg_hh.get_lost().items():
    print("\t- {} at pg/nv/hh has been lost in ".format(hog))
    if len(genomes) == 1:
    	for g in genomes:
    		if format("{}".format(g)) == "Halyomorpha_halys":
    			hh_hogs_lost.append(hog)
print("\n")

hh_lost = open("hh_vs_pg_lost", "w")
for h in hh_hogs_lost:
	lost_genes = h.get_all_descendant_genes()
	for element in lost_genes:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		lost_dict = gene_name.get_dict_xref()
		lost_genes = lost_dict['protId']
		hh_lost.write(lost_genes+"\n")	


#now get the genes lost in Euschistus heros
eh_hogs_lost = []
print("LOST GENES")
for hog, genomes in pg_eh.get_lost().items():
    print("\t- {} at pg/nv/hh has been lost in ".format(hog))
    if len(genomes) == 1:
    	for g in genomes:
    		if format("{}".format(g)) == "Euschistus_heros":
    			eh_hogs_lost.append(hog)
print("\n")

eh_lost = open("eh_vs_pg_lost", "w")
for h in eh_hogs_lost:
	lost_genes = h.get_all_descendant_genes()
	for element in lost_genes:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(element))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		lost_dict = gene_name.get_dict_xref()
		lost_genes = lost_dict['protId']
		eh_lost.write(lost_genes+"\n")
		



#Make ancestral dataframe

def get_ancestral_hog_id(ancestral_hog):
    '''Ancestral hog id for this is a concatenation of all descendant genes'''
    ancestral_hog_id = str(ancestral_hog.get_all_descendant_genes())
    return ancestral_hog_id
    	
###get ancestral genome

def make_ancestral_dataframe(ancestral_genome_name, hamObj):
    '''
    This is where all the information about a particular ancestral genome is contained.
    I name the ancestral genes as the concatenation of all the descendant genes.
    Each row is a different descendant gene, thus there may be multiple rows with the same ancestral gene.
    '''
    #get ancestral genome using pyham
    ancestral_genome = hamObj.get_ancestral_genome_by_name(ancestral_genome_name)

    #get all ancestral genes (hogs)
    ancestral_genes = ancestral_genome.genes
    ancestral_genome_df = pd.DataFrame(ancestral_genes, columns={"ancestral_gene"})

    #get descendant genes of the hog
    ancestral_genome_df['descendant_genes'] = ancestral_genome_df.apply(lambda x: x["ancestral_gene"].                                                                        get_all_descendant_genes(), axis=1)

    #modify df so that each descendant gene has its own row
    ancestral_genome_df = ancestral_genome_df.set_index(['ancestral_gene'])['descendant_genes'].apply(pd.Series).stack()                                             .reset_index(level=1, drop=True).reset_index()

    #rename the column
    ancestral_genome_df = ancestral_genome_df.rename({0: "descendant_gene"}, axis=1)

    #make a string id of ancestral gene
    ancestral_genome_df['ancestral_gene_id'] = ancestral_genome_df.apply(lambda x: get_ancestral_hog_id(x['ancestral_gene']),                                                                         axis=1)

    #convert pyham gene name to cross-reference gene name
    ancestral_genome_df['descendant_gene_xref'] = ancestral_genome_df.apply(lambda x: x['descendant_gene'].                                                                            get_dict_xref()['protId'].                                                                            split(" ")[0], axis=1)
    
    print("There are {} ancestral genes in the {} genome.".          format(len(ancestral_genome_df.groupby('ancestral_gene_id').size()), ancestral_genome_name))

    return ancestral_genome_df

an_df = make_ancestral_dataframe("Piezodorus_guildinii/Nezara_viridula/Halyomorpha_halys", ham_analysis)


an_df.to_csv("ancestral_df.csv", sep="\t")

'''ancestral_genes = pg_nv_hh.get_ancestral_clustering()
#print(ancestral_genes)
print(len(ancestral_genes))

an_genome = open("ancestral_genes.txt" ,"w")
for h,v in ancestral_genes.items():
	for g in v:
		gene_num = re.findall(r"\(\s*\+?(-?\d+)\s*\)", str(g))
		gene_name = ham_analysis.get_gene_by_id(gene_num[0])
		an_dict = gene_name.get_dict_xref()
		an_genes = an_dict['protId']
		an_genome.write(an_genes+"\n")

#print(ancestral_genes)'''	
	
	
lost_of.close()
gained_of.close()
dup_of.close()
nv_lost.close()
hh_lost.close()
eh_lost.close()
#an_genome.close()




