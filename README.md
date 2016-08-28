
# yMap - Yeast Genotype to Phenotype Map (release year 2016)

yMap is a python based fast and robust automated method to map large yeast variants data to 
            
            - proteins post-translational modifications 
                            
            - proteins domains, 
                            
            - proteins-nucleotide binding domains, 
                        
            - proteins structural regions, 
                
            - proteins active and binding sites 
                
            - proteins networks visualisation. 


The post-translational modifications in yMap are collected from different repositories like UniProt and sources 
with annotated PTMs like PTMcode 2.0 and PTMfunc, for more details, see below.

In a user friendly three steps, it generates a "final-report" file to report all the non-synonymous 
mutations that overlaps or falls inside the above mentioned proteins functional regions.
The final-report is complemented with two other files; enrichment and visualsation id file

#Dependencies

yMap depends on:

        python 2.6.x

        python 3.x

        Orange bioinformatics (http://pythonhosted.org/Orange-Bioinformatics/#installation)

#Installation

            pip install ymap

#Usage
    
            step1: $ ydata       	#download all the data needed for proper execution of ymap
            
            step2: copy and paste the "mutation file" to the present directory

            step3: $ yproteins  	#if starting file contains the mutations at proteins level 
                                        (SEE example_mutation_file/mutation.txt).
            
            step3: $ ygenes        	#if starting file contains the mutations at chromosomes leve with genetic                                                    coordinates (SEE example_mutation_file/mutated_proteins.txt).
            
            step4: $ yweb		 # generates the html based visualization of mutated proteins on BioGrid db.
					(NOTE: a user will required to specify the 'path/to/biog.txt' as input, when asked)
					

*To run from source code:

		Change path to directory containing ymap.py
		
		$python ymap.py -d ydata (step1)
		
		$pyhton ymap.py -p yproteins (step3)
		
		$pyhton ymap.py -g ygenes (step3)
		
		$pyhton ymap.py -w yweb (step4)

#Contents:
Introduction to different types of data (generated/provided in yMap)
Introduction to all the methods
Results (introduction to results data)
Troubleshoots

#Introduction to data:

	—————input———

A - mutation (tab separated txt "mutated_proteins.txt") file contains proteins common names and mutated residues positions
(please following the exact naming convention of input files as in example data, for proper execution of ymap; see example data))


	———output———(Pre-analysis data needed for ymap execution)


(i)	Raw files downloaded from UniProt and stored in the present dir.by executing step2.  	


1 - uniprot_mod_raw.txt	# Uniprot data in raw format 


2 - yeastID.txt			# Yeast id containing file


3 - PTMs.txt			# contains yeast proteins, PTMs position and PTM types


4 - PTM_id_file.txt		# combined file of 2 and 3.


5 - domains.txt			# yeast proteins, domains start, end and names


6 - id_domain.txt		# combined file of 2 and 5.


7 - bact.txt			# contains proteins id, and binding and active sites 
				  positions

8 - sites_id.txt		# combined file of 2 and 7.


9 - uniprot_bioGrid.txt		# contains all the yeast proteins with BioGrid ids


(i-B)	Pre downloaded files from PTMcode and PTMfunc

(PTMfunc)

3DID_aceksites_interfaceRes_sc.txt

3DID_phosphosites_interfaceRes_sc.txt

3DID_ubisites_interfaceRessc_sc.txt

SC_psites_interactions_sc.txt

SC_ubi_interactions_sc.txt

SC_acet_interactions.txt

(PTMcode)

schotspot.txt

sc_btw_proteins.txt

sc_within_proteins.txt



(ii)	Processed data from UniProt and other resources by executing step2.  

A number of files germinated from the original UniProt file for further analyses:

PTMs.txt			
		contains Post-translational modifications

PTM_id_file.txt			
		PTMs.txt with all the proteins ids

PDB.txt				
		contains PDB structural data from UniProt

nucleotide.txt			
		contains DNA-Protein binding motifs 

back.txt			
		contains Proteins active and binding positions

d_id_map.txt			
		contains protein domains with all the ids	

id_domain.txt			
		gff data from frmt.txt with all the ids

domains.txt			
		domains data from UniProt

frmt.txt			
		formatted gff file for further process

sites_id.xt			
		Active/binding sites with all ids

unipro_bioGrid.txt		
		contains BioGrid ids of all yeast proteins

nucleotide.txt			
		proteins (uniprot) id with DNA binding motifs

id_nucleotide.txt		
		contains data from nucleotide with all the protein ids for processing

#Results

(inside ymap-results folder, each subfolder contains three files, one with mutations analysis file, which includes mutated proteins, mutation positions, mutated functional region and source of data, pvalue.txt of pathways enrichments and biog.txt, a biogrid id corresponding to mutated proteins)

	/PTMs/mutated_proteins.txt	
		contains proteins ids mutated at PTMs sites

	/Domains/domains_mapped.txt	
		contains proteins ids mutated for protein domains

	/A-B_binding/ab_mutation_file.txt	
		contains proteins ids mutated at active and binding 


	PPI - PTMfunc data

		PPI/acetylation
		PTM-type containing residue is important in PPI

		PPI/Phosphorylation
		PTM-type containing residue is important in PPI

		PPI/ubiquitination
		PTM-type containing residue is important in PPI

	Interface

		Interface/ubiquitination
		PTM-type containing residue present at protein interface 

		Interface/acetylation
		PTM-type containing residue present at protein interface 

		Interface/Phosphorylation
		PTM-type containing residue present at protein interface 

	PTMs_hotSpot
		PTMs concentrated in a small motif known as hopspot by Beltrao et al. Cell 2012.

	PTMs_between_proteins - PTMcode2.0 data
		PTMs present between two proteins and involvined in crosstalk. 

	PTMs_witnin_proteins
		PTMs present within a protein and involvined in crosstalk.
	
	biog.txt			
		contains proteins BioGrid ids for -w web function (this file present in each subfolder).

	p-value.txt			
		contains pathways enrichments for each type of mutation observed (this file present in each subfolder).

	summary.txt			
		contains all the proteins that are mapped on different data sets.

	final_report.txt		
		its a refined version of summary.txt and contains, protein UniProt id, common names, amino acid mutation position, wild type amino acid, mutated amino acid, type of mutation (non-synonymous/stop codon), mutation feature types (i.e. PTM-type or domain-name etc), mutation feature (i.e. PTMs, domain or another) and source of data (e.g. UnProt)



#Introduction to all the methods 
		(How individual methods work in ymap)

NOTE: change the name of the mutations containing file to ‘mutated_proteins.txt’ (see example data) and copy to the cd path/to/ymap


Functions name			Description

mutation_types_file()	mutation type and amino acid change calculation (where ref. and mutant base known)

pTMdata()		
		Downloads UpiProt data as a raw txt file (uniprot_mod_raw.txt)

clean()			
		cleans file 'uniprot_mod_raw.txt' into a tab separated’PTMs.txt'

iD()			
		This method retrieves the different ID types for maping (yeastID.txt)

pmap()			
		if proteins ids are not SDG or uniprot or common names, this method maps the ids 

ptm_map()		
		This method maps the overlap between mutated codons from previous method to the PTM sites 

dclean()		
		domain data needed to be filters from UniProt file, before mapping domains 
d_map()			
		maps mutations to the yeast domains (id_domain.txt)

dmap()			
		map mutations to proteins domains (domains_mapped.txt)

enrich()		
		This method performed enrichment analysis of mutated proteins and return the p value of functional enrichment 
		of mutated proteins at different functional regions/residues; see main text for how pvalue is calculated. 
ab()			
		prepares raw Uniprot data (uniprot_mod_raw.txt) for yeast active and binding sites mutation analysis (bact.txt)

id()			
		maps proteins ids to active and binding sites containing proteins (sites_id.txt)

mmap()			
		map mutations to proteins active and bindings sites (ab_mutation_file.txt)

nucleotide()		
		prepares the UniProt data for the nucleotide motifs mapping to mutations

n_map()		
		maps different proteins ids to nucleotides data

nucleotide_map()	
		maps mutations to the nucleotide binding motifs

bioGrid()		
		Downloads BioGrid ids of yeast proteins from UniProt for further processing including mapping and web browsing
        WARNING: requires powerful machines to work with as its expensive to open in machines with low memory.

preWeb()
		maps mutations to BioGrid ids (biog.txt)

bweb()		opens the BioGrid db in browser with as many tabs as mutated proteins

pdb_c()		Structure data filtration from UniProt

mu_map()	mutations proteins mapped to the yeastID file

pdb()		This code maps mutations to the proteins structural regions

interface()	PTM present at the interface of two proteins and known to play role in interaction (Beltrao et al. Cell 2012)

ppi()		PTM present at the interface of two proteins and known to play role in interaction (Beltrao et al. Cell 2012)

withinPro()	PTMs (predicted) involved in the crosstalk within a given protein at baker's years (Minguez el 2012)

betweenPro()	PTMs (predicted) involved in the crosstalk in different proteins at baker's years (PTMcode 2.0; Minguez el 2012)

hotspot()	PTMs containing motifs in a close proximity are named hotspots (Beltrao et al. Cell 2012)

#Troubleshoots


1 - The files of annotated PTMs are missing or less them nine.

Reason: unzip the data/PTMcode+PTMfunc_data/sc_btw_proteins.txt.zip did not worked in $ ydata command.
how to correct: manually unzip the sc_btw_proteins.txt.zip file and run $ ydata (normally this will not needed)

2 - $ ygenes gives an error message: 

“IndexError: string index out of range”

2(b) - The same reason (below) leads to the unsuccessful mapping of mutations to different functional regions like domains:

"Error: input file contains error position forBRR2protein"

Reason: the mutations positions fall outside the start and end of the respective proteins (NOTE: to analyse
the proteins in starting file with correct mutation positions, user can use individual methods uniprot_data()
and functional_data(), to get all the analyses done, than execute the command-line step3)

how to correct: Look at the positions of mutations and compare them manually if they correspond to start and end 
positions of a protein, if not, correct the problem and re-run $ ygenes command.

3 - yweb fails to locate the directory.

how to correct: In python 2.x, the path should be given as “path/to/biog.txt” but in python 3.x it’s without inverted commas, 
path/to/biog.txt


# Contributors

            http://www.biw.kuleuven.be/CSB/
            
This work in supported by KU Leuven research fund. 
