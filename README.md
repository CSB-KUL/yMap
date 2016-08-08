
# yMap - Yeast Genotype to Phenotype Map (release year 2016)

yMap is a python based fast and robust automated method to map large yeast variants data to 
            
            - proteins post-translational modifications 
                            
            - proteins domains, 
                            
            - proteins-nucleotide binding domains, 
                        
            - proteins structural regions, 
                
            - proteins active and binding sites 
                
            - proteins networks visualisation. 


The post-translational modifications in yMap are collected from different repositories like UniProt and sources 
with annotated PTMs like PTMcode 2.0 and PTMfunc, for more details, see README file.

In a user friendly way, it generates a "final-report" file to report all the non-synonymous 
mutations that overlaps or falls inside the above mentioned proteins functional regions.
The final-report is complemented with two other files; enrichment and visualsation id file

#Dependencies

yMap depends on:

        python 2.7.x

        python 3.5.x

        Orange bioinformatics (http://pythonhosted.org/Orange-Bioinformatics/#installation)

#Installation

            pip install ymap

#Usage
    
            cd /path/to/ymap
            
            $ python ymap.py --data data        #download all the data need for proper execution of ymap
            
            $ python ymap.py --ymap_proteins ymap_proteins  #if starting file contains the mutations at proteins level 
                                                            (SEE example data).
            
            $ python ymap.py --ymap_genes ymap_genes        #if starting file contains the mutations at chromosomes leve with genetic                                                    coordinates (SEE example data).
            
            $ python ymap.py --web web --web_path path   # generates the html based visualization of mutated proteins on BioGrid db.
		(NOTE: a user will required to specify the path/to/biog.txt as input, when asked)

#Contents:
Introduction to different types of data (generated in yMap)
Introduction to all the methods
Manual (steps of using yMap)
Troubleshoots

#Introduction to data:

	—————input———

A - mutation (tab separated txt (test)”mutated_proteins.txt”) file contains proteins common names and mutated residues positions(see example data))


	———output———


(i)	Raw files (downloaded from UniProt and stored in the present dir.) 	


1 - uniprot_mod_raw.txt	# Uniprot data in raw format 


2 - yeastID.txt			# Yeast id containing file


3 - PTMs.txt			# contains yeast proteins, PTMs position and PTM types


4 - PTM_id_file.txt		# combined file of 2 and 3.


5 - domains.txt			# yeast proteins, domains start, end and names


6 - id_domain.txt		# combined file of 2 and 5.


7 - bact.txt			# contains proteins id, and binding and active sites 
				  positions

8 - sites_id.txt		# combined file of 2 and 7.


9 - uniprot_bioGrid.txt	# contains all the yeast proteins is and BioGrid ids


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



(ii)	Processed data from UniProt and other resources. 

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

(each folder contains three files, one with mutations analysis file, which includes mutated proteins, mutation positions, mutated functional region and source of data, pvalue.txt of pathways enrichments and biog.txt, a biogrid id corresponding to mutated proteins)

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
		contains proteins BioGrid ids for web() function

	p-value.txt			
		contains pathways enrichments for each type of mutation observed. 

	summary.txt			
		contains all the proteins that are mapped on different data sets.

	final_report.txt		
		contains, protein UniProt id, common names, amino acid mutation position, wild type amino acid, mutated 
		amino acid, type of mutation (synonymous/ non-synonymous/stop codon), mutation feature types (i.e. PTM-type or domain-name etc), mutation feature (i.e. PTMs, domain or another) and source of data (e.g. UnProt)



#Introduction to all the methods 
		(How individual methods work in ymap)

NOTE: change the name of the mutations containing file to ‘mutated_proteins.txt’ (see example data) and copy to the cd path/to/ymap


Functions name			Description

mutation_types_file()		
		to get the mutations quality (Synonymous, Non-Synonymous or Stop codon) in the file containing genome coordinates of mutation with knowledge of ref. and mutated base; furthermore the ‘mutation.txt’ from this code will be used to further other methods n map to get the summary file.

pTMdata()		
		Downloads UpiProt data as a raw txt file (uniprot_mod_raw.txt)

clean()			
		cleans file 'uniprot_mod_raw.txt' into 'PTMs.txt'

iD()			
		Downloads different yeast IDs from UniProt (yeastID.txt)

pmap()			
		combine yeast PTM and IDs into one file(PTM_id_file.txt)

ptm_map()		
		Gives a file containing mutated proteins at PTMs with their positions 

dclean()		
		prepares raw Uniprot data (uniprot_mod_raw.txt) for yeast domains mutations analysis (domains.txt)

d_map()			
		combine domains.txt and yeastID.txt into one file (id_domain.txt)

dmap()			
		map mutations to proteins domains (domains_mapped.txt)

enrich()		
		Gives p-value of pathways underlying the mutated proteins at for all type of data analysed the methods (either
		PTMs, domains and/or active - binding sites, structural and DNA-protein binding motifs or PTMcode/PTMfunc data)

ab()			
		prepares raw Uniprot data (uniprot_mod_raw.txt) for yeast active and binding sites mutation analysis (bact.txt)

id()			
		combine bact.txt and yeastID.txt into one file (sites_id.txt)

mmap()			
		map mutations to proteins active and bindings sites (ab_mutation_file.txt)

nucleotide()		
		prepares the UniProt data for the nucleotide motifs mapping to mutations

n_map()		
		maps the file from previous method to yeastID file

nucleotide_map()	
		maps mutations to the nucleotide binding motifs

bioGrid()		
		Downloads BioGrid ids of yeast proteins from UniProt (uniprot_bioGrid.txt)

preWeb()
		combines mutated proteins with uniprot_bioGrid.txt in a file (biog.txt)

web()			Automated opening the web links of BioGrid db for further network and pathways analyses of mutated proteins. WARNING. It requires powerful machines to open new tabs for each mutated protein.

ymap()			runs all the methods included in ymap (excluding data() and mutation_types_file())  and returns different data folders for each type of data processed for our mutational data.



#Troubleshoots


1 - The files of annotated PTMs are missing or less them nine.

Reason: forget to unzip the data/PTMcode+PTMfunc_data/sc_btw_proteins.txt.zip
how to correct: unzip the sc_btw_proteins.txt.zip file and run the data()

2 - $ python ymap.py --ymap_genes genes gives an error message: 

“IndexError: string index out of range”

Reason: the mutations positions fall outside the start and end of the respective proteins (NOTE: to analyse
the proteins in starting file with correct mutation positions, user can use individual methods uniprot_data()
and functional_data(), to get all the analyses done)

how to correct: Look at the positions of mutations and compare them manually if they correspond to start and end 
positions of a protein, if not, correct the problem and re-run the method mutation_types_file().


# Contributors

            http://www.biw.kuleuven.be/CSB/
            
This work in supported by KU Leuven research fund. 
