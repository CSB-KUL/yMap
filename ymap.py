#!/usr/bin/python
#
# year of release ! 2016
#
# With the exponential growth of Posttranslational modifications (PTMs) 
# data and and lack of characterisation of all the PTM-types. Its important
# to undersand their functions and experimental relevence. 
# And to understand the importance of PTMs in Yeast genome, its important 
# to make it easier to map mutations to PTM positional data. Besides that
# it also important to translate genetic abrretion to understand the phenotype. 
# We architect the yMap library to help users to understand which parts of
# mutated proteins are affected during the yeast experimentation.
# This facilitation not only would help bioligists to interpret their data 
# efficiently but also gives freedom to save time by mapping data to mutations 
# easily 
#
# The yMap program is a python based fast and robust automated method to map 
# large yeast variants to proteins post-translational modifications, proteins domains,
# proteins-DNA binding domains, proteins structural regions, proteins active and 
# binding sites, proteins networks visualisation. 
#
# Dependencies:
# Orange Bioinformatics 
# pip install orange; see README file 
#
# For Usage see README file

import os
from itertools import groupby
import shutil
import time
import urllib
import urllib2
from collections import OrderedDict
import webbrowser
import orange
from orangecontrib.bio import go    
ontology = go.Ontology()
annotations = go.Annotations("sgd", ontology=ontology)  


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#/////////////////////// Codon calculation in case mutation position is not in proteins ///////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class codon_from_protein_file:

    def gff(self):
        """ The genomic coordinates downloaded in gff formate for further processing to calculate mutated codons, if not
              available, see next method"""
        url = 'http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff'
        request = urllib2.Request(url)
        request.add_header('User-Agent', 'Python %s')
        response = urllib2.urlopen(request)
        page = response.read()
        file = open('gff.txt','w')  
        file.write(page)
        file.close()
        

    def frmt(self, file_gff):
        """This method format the gff file into a tsv one, with protein id, start and end with strand orientation"""
        with open('frmt.txt','w') as file4:
            with open(file_gff, 'r') as file_gff:
                for line in file_gff:
                    if not line.startswith('##') and not line.startswith('#'):
                        word = line.split()
                        if len(word)!=1 and word[2]=='gene':
                            result = word[3]+'\t'+word[4]+'\t'+word[6]+'\t'+word[8]
                            result = result.split()
                            result = result[3].split(';')
                            results = result[0].split('=')
                            result2 = results[1]+'\t'+word[3]+'\t'+word[4]+'\t'+word[6]
                            file2 = open('frmt.txt','a')
                            file2.write(result2+'\n')
                    
    
    def ID(self):

        """ This method retrieves the different ID types for maping """

        url='http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,genes(OLN),%2Cgenes(PREFERRED)'
        request = urllib2.Request(url)
        request.add_header('User-Agent', 'Python %s')
        response = urllib2.urlopen(request)
        page1 = response.read()
        file =open('yeastID.txt','w')       
        file.write(page1)
        file.close()

    def id_map(self, file_id_name, frmt):        
        with open('d_id_map.txt', 'w') as file2:
            with open(file_id, 'r') as file_id_name:
                for line in file:
                    line=line.split()
                    with open(frmt, 'r') as fr:
                        for f in fr:
                            f=f.split()
                            if len(line)>2:
                                if line[1]==f[0]:
                                result= line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+f[1]+'\t'+f[2]+'\t'+f[3]
                                if result>0:
                                    file2= open('d_id_map.txt', 'a')
                                    file2.write(result+'\n')



    def cc_mapping(self, mu, frmt_file):
        """continuation of codon calculation (cc) with formated gff file, the file product is the protein id with mutated codons"""
       with open('mutation.txt', 'w') as file4:
            with open(frmt_file,'r') as f_file:
                for line in f_file:
                    line=line.split()
                    with open(mu, 'r') as m:
                        for m in m:
                            m=m.split()
                            if line[2]==m[0]:
                                q= int(float(line[4]))-int(float(m[1]))
                                q1=int(q)/3
                                q2=m[0]+'\t'+str(q1)
                                file4 = open('mutation.txt', 'a')
                                file4.write(q2+'\n')
                                

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Mutation type (Synon | Non-Synon | Stop codon) module (see exmple data) \\\\\\\\\\\\\\\\\\\\\
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


def translate_dna(dna):
    """ calculate the start position for the final codon """
    last_codon_start = len(dna) - 2 
    protein = "" 
    # process the dna sequence in three base chunks
    for start in range(0,last_codon_start,3): 
        codon = dna[start:start+3]
        aa = genetic_code.get(codon.upper(), 'X') 
        protein = protein + aa 
    return protein 


def revcomp(dna, reverse=True, complement=True):
    """ reverse complement of a protein in negative strand"""
    bases = 'ATGCTACG'
    complement_dict = {bases[i]:bases[i+4] for i in range(4)}
    if reverse:
        dna = reversed(dna)
        result_as_list = None
    if complement:
        result_as_list = [complement_dict[base] for base in dna]
    else:
        result_as_list = [base for base in dna]
    return ''.join(result_as_list)


def mutation_file(file1, file2):
    """  The following bit of code combines, translation and reverse complement codes and out put the 
         type of mutation, amino acid, both wild and mutant. 
    """
    with open('mutation.txt', 'w') as t:
        with open(file1, 'rU') as mut:
            for m in mut:
                m = m.split()
                with open(file2,'r') as id:
                    for i in id:
                        i = i.split()
                        if not m[0].startswith('c'.upper()):
                            if len(m) != 5  or not m[0].startswith('c'.lower()):
                                raise StopIteration('Please correct the format of input mutation file')
                            else:
                                if m[4] == i[2]:
                                    take = m[4]+'\t'+m[0]+'\t'+i[3]+'\t'+m[1]+'\t'+i[4]+'\t'+m[2]+'\t'+m[3]+'\t'+i[5]
                                    take1= take.split()
                                    cod = int(take1[4])-int(take1[3])
                                    c = int(cod)/3
                                    with open('gff.txt', 'rU') as orf:  
                                        linee = orf.readlines()[23078:]
                                        up = (x[1] for x in groupby(linee, lambda line: line[0] == ">"))
                                        for head in up:
                                            head = head.next()[1:].strip()
                                            seq = "".join(s.strip() for s in up.next())
                                            if head == take1[1] and take1[0] == i[2] and take1[7] == '-':
                                                cn = c-1
                                                sli_n = seq[int(take1[2]):int(take1[4])]
                                                sli_m_n = seq[int(take1[2]):int(take1[3])]+take1[6]+ seq[int(take1[3])+1:int(take1[4])]
                                                rev_sli_n = revcomp(sli_n, reverse=True, complement=True)
                                                rev_sli_m_n = revcomp(sli_m_n, reverse=True, complement=True)
                                                wild_type_rev_n = translate_dna(rev_sli_n)
                                                mut_type_n = translate_dna(rev_sli_m_n)
                                                if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] == '_':
                                                    pic = take1[0]+'\t'+str(c)+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Stop'
                                                    if pic > 0:
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic+'\n')
                                                if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] != '_':
                                                    pic = take1[0]+'\t'+str(c)+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Non-Synonymous'
                                                    if pic > 0:
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic+'\n')
                                                        continue
                                            if head == take1[1] and take1[0]==i[2] and take1[7] == '+':
                                                cp = c-1
                                                pos = int(take1[2])-1
                                                sli_p = seq[int(pos):int(take1[4])]
                                                sli_m_n = seq[int(pos):int(take1[3])]+take1[6]+seq[int(take1[3])+1:int(take1[4])]
                                                wild_type_p = translate_dna(sli_p)
                                                mut_type_p = translate_dna(sli_m_n)
                                                if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp] != '_':
                                                    pick = take1[0]+'\t'+str(c)+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Non-Synonymous'
                                                    if pick > 0:
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pick+'\n')
                                                if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp]=='_':
                                                    pick = take1[0]+'\t'+str(c)+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Stop'
                                                    if pick > 0:
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pick+'\n')


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# //////////////////                     UniProt data                 /////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class YGtPM:

    def PTMdata(self):

        """download modified residues data from uniprot and save in tsv frmt text time"""

        url = 'http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=gff&columns=id,feature(MODIFIED%20RESIDUE)'
        request = urllib2.Request(url)
        request.add_header('User-Agent', 'Python %s')
        response = urllib2.urlopen(request)
        page = response.read()
        file = open('uniprot_mod_raw.txt','w')  
        file.write(page)
        file.close()
        return "PTM raw file is ready for formting" 
    
    def clean(self, UniProt_file):              
        
        """ This method clean the downloaded PTMs data from web into a tsv file
            UniProt_file_name = uniprot_mod_raw.txt
        """

        with open('PTMs.txt', 'w') as out:
            with open(UniProt_file,'r') as UniProt_file_name:
                for l in UniProt_file_name:
                    if not l.startswith('##'):
                        line = l.split()
                        if line[2] == 'Lipidation':
                            lll = line[0]+'\t'+line[4]+'\t'+line[8]
                            ll = lll.split()
                            ll = ll[2].split('=')
                            p = line[0]+'\t'+line[4]+'\t'+ll[1]
                            out.write(p+'\n')
                            continue
                        if line[2] == 'Glycosylation':
                            ggg = line[0]+'\t'+line[4]+'\t'+line[8]
                            gg = ggg.split()
                            gg =  gg[2].split('=')
                            pp = line[0]+'\t'+line[4]+'\t'+gg[1]
                            if pp > 0:
                                out = open('PTMs.txt', 'a+')
                                out.write(pp+'\n')
                                continue
                        if line[2] == 'Modified':
                            mmm = line[0]+'\t'+line[4]+'\t'+line[9]
                            mm = mmm.split()
                            mm = mm[2].split('=')
                            mm = mm[1].split(';')
                            ppp = line[0]+'\t'+line[4]+'\t'+mm[0]
                            if ppp > 0:
                                out = open('PTMs.txt', 'a+')
                                out.write(ppp+'\n')
                                continue
                        if line[2] == 'Cross-link': #ubiquitination
                            ccc = line[0]+'\t'+line[4]+'\t'+line[8]
                            cc = ccc.split()
                            cc = cc[2].split('=')
                            pppp = line[0]+'\t'+line[4]+'\t'+cc[1]
                            if pppp > 0:
                                out = open('PTMs.txt', 'a+')
                                out.write(pppp+'\n')


    def ID(self):

        """ This maethod retrieves the different ID types for maping """

        url='http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,genes(OLN),%2Cgenes(PREFERRED)'
        request = urllib2.Request(url)
        request.add_header('User-Agent', 'Python %s')
        response = urllib2.urlopen(request)
        page1 = response.read()
        file_1 =open('yeastID.txt','w')  
        file_1.write(page1)
        file_1.close()
        
            
    def Map(self, file_id, file_PTMs):          #file_PTMs = PTMs.txt from clean()

        """ if proteins ids are not SDG or uniprot or common names, this method maps the ids"""

        with open('PTM_id_file.txt', 'w') as file3:
            with open(file_id, 'r') as file_id_name:
                for lin in file_id_name:
                    line = lin.split()
                    with open(file_PTMs) as ptms:
                        for i in ptms:
                            i = i.split()
                            if len(line) > 2:
                                if line[0] == i[0]:
                                    result3 = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+i[1]+'\t'+i[2]
                                    if result3 > 0:
                                        file3 = open('PTM_id_file.txt', 'a')
                                        file3.write(result3+'\n')
                                             

    def ptm_map(self, mutation_file, PTM_id_file):

        """ This method maps the overlap between mutated codons from previous method to the PTM sits"""
        summary = open('summary.txt', 'w')
        with open('mutated_proteins.txt', 'w') as file5:
            with open(mutation_file, 'r') as mutation_file:
                for line in mutation_file:
                    line = line.split()
                    with open(PTM_id_file, 'r') as file_PTMs:
                        for line1 in file_PTMs:
                            line1 = line1.split()
                            if line[0] == line1[2] and line[1] == line1[3]:
                                take = line1[0]+'\t'+line1[1]+'\t'+line[0]+'\t'+line[1]+'\t'+line1[4]+'\t'+'UniProt'
                                if take > 0:
                                    file5 = open('mutated_proteins.txt', 'a') # contains proteins, PTMs mutation site, PTM-types in txt file
                                    summary = open('summary.txt', 'a')
                                    file5.write(take+'\n')
                                    summary.write(line1[0]+'\t'+line[0]+'\t'+line[1]+'\t'+line1[4]+'\t'+'PTMs'+'\t'+'UniProt'+'\n')


    def dclean(self, domain):  

        """domain data needed to be filters from UniProt file, before mapping domain = uniprot_mod_raw.txt from 
        previous method"""

        with open('domains.txt', 'w') as file_d:
            with open(domain, 'r') as d:
                for line in d:
                    if not line.startswith('##'):
                        line = line.split()
                        if line[2] == 'Domain':
                            take = line[0]+'\t'+line[3]+'\t'+line[4]+'\t'+line[8]
                            take1 = take.split()
                            take1 = take1[3].split(';')
                            take2 = line[0]+'\t'+line[3]+'\t'+line[4]+'\t'+take1[0]
                            if take2 > 0:
                                file_d = open('domains.txt', 'a')
                                file_d.write(take2+'\n')
                
    
    def d_map(self, yeast_id, domain):

        """ For multi proteins Ids, this code maps proteins UniProt ids to multiple Ids
            yeast_id=yeastID.txt """  

        id_domain = open('id_domain.txt', 'w')
        with open(yeast_id, 'r') as fl:
            for f in fl:
                f = f.split()
                with open(domain,'r') as dp:
                    for d in dp:
                        d = d.split()
                        if len(f) > 2:
                            if f[0] == d[0]:
                                take =d [0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]
                                if take > 0:
                                    id_domain=open('id_domain.txt', 'a')
                                    id_domain.write(take+'\n')


    def dmap(self, file1, file2): 

        """ maps mutations to domains regions #file2 = id_domain.txt """  

        with open('domains_mapped.txt', 'w') as mp:
            with open(file1,'r') as f:
                for line in f:
                    line1=line.split()
                    with open(file2, 'r') as f2:
                        for line2 in f2:
                            line2 = line2.split()
                            if line1[0] == line2[2]:
                                if int(line1[1]) >= int(line2[3]) and int(line1[1]) <= int(line2[4]):
                                    take = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+'UniProt'
                                    if take > 0:
                                        mp= open('domains_mapped.txt', 'a')
                                        summary = open('summary.txt', 'a+')
                                        mp.write(take+'\n')
                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+'domains'+'\t'+'UniProt'+'\n')

 
    def enrich(self, file1):

        """ This method performed enrichment analysis of mutated proteins and
        return the p value of functional enrichment of mutated proteins at PTMs residues; 
        file1 = mutated_proteins.txt OR domains_mapped.txt"""

        with open('pvalue.txt','w') as out:
            with open(file1, 'r') as f:
                for line in f:
                    line=line.split()
                    res = annotations.get_enriched_terms(line)
                    if len(res) == 0:
                        out = open('pvalue.txt','a')
                        out.write('No enrichment found')
                    else:
                        for go_id, (genes, p_value, ref) in res.items():
                            if p_value < 0.002:
                                take = ontology[go_id].name, p_value, ''.join(genes)
                                out.write(str(take))
                                continue
                            if p_value < 0.005:
                                take1 = ontology[go_id].name, p_value, ''.join(genes)
                                out.write(str(take1))
                                continue
                            if p_value < 0.02:
                                take2 = ontology[go_id].name, p_value, ''.join(genes)
                                out = open('pvalue.txt','a+')
                                out.write(str(take2))
                                continue
                            if p_value < 0.05:
                                take3 = ontology[go_id].name, p_value, ''.join(genes)
                                out = open('pvalue.txt','a+')
                                out.write(str(take3))


    def ab(self, file_raw): 

        """Proteins active and binding sites mapping (uniprot_mod_raw.txt) from PTMdata() 
        ab() cleans the data; file_raw='uniprot_mod_raw.txt' """

        with open('bact.txt','w') as file2:
            with open(file_raw, 'r') as d:
                for f in d:
                    if not f.startswith('##'):
                        f = f.split()
                        if f[2] == 'Active':
                            take = f[0]+'\t'+f[2]+'\t'+f[4]
                            if take > 0:
                                file2 = open('bact.txt','a')
                                file2.write(take+'\n')
                        if f[2] == 'Binding':
                            take2 = f[0]+'\t'+f[2]+'\t'+f[4]
                            if take2 > 0:
                                file2 = open('bact.txt','a+')
                                file2.write(take2+'\n')
                            
                             
    def id(self, act, yeast_id): 

        """ perform only if you have mutation list with common gene names 
        #act = bact.txt; yeast_id=yeastID.txt
        """
        with open('sites_id.txt', 'w') as file_id:
            with open(act, 'r') as a:
                for a in a:
                    a = a.split()
                    with open('yeastID.txt')as id:
                        for i in id:
                            i = i.split()
                            if len(i) > 2:
                                if a[0] == i[0]:
                                    take = i[2]+'\t'+i[1]+'\t'+i[0]+'\t'+a[1]+'\t'+a[2]
                                    if take > 0:
                                        file_id = open('sites_id.txt', 'a')
                                        file_id.write(take+'\n')


    def mmap(self, file_sites, mutation):

        """ maps mutations to ab (active and binding sites) file_sites= sites_id.txt""" 

        with open('ab_mutation_file.txt', 'w') as out: 
            with open(file_sites, 'r') as s:
                for a in s:
                    a = a.split()
                    with open(mutation, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if a[0] == m[0] and a[4] == m[1]:
                                take = a[2]+'\t'+ a[3]+'\t'+m[1]+'\t'+'UniProt'
                                if take > 0:
                                    out = open('ab_mutation_file.txt', 'a')
                                    summary = open('summary.txt', 'a+')
                                    out.write(take+'\n')
                                    summary.write(a[2]+'\t'+a[0]+'\t'+m[1]+'\t'+ a[3]+'\t'+'Active/Binding site'+'\t'+'UniProt'+'\n')


    def nucleotide(self):

        """ Proteins nucleotide binding motifs mapping (uniprot_mod_raw.txt) from PTMdata()
        data filteration from UniProt file for DNA-Proteins motifs"""

        with open('nucleotide.txt', 'w') as t:
            with open('uniprot_mod_raw.txt', 'rU') as file_raw:
                for fi in file_raw:
                    if not fi.startswith('##'):
                        f = fi.split()
                        if f[2] == 'Nucleotide' and len(f) > 8:
                            take = f[0]+'\t'+f[2]+'\t'+f[4]+'\t'+f[5]+'\t'+f[9]
                            take1 = take.split()
                            take1 = take1[4].split(';')
                            if take > 0:
                                t = open('nucleotide.txt', 'a')
                                t.write(f[0]+'\t'+f[2]+'\t'+f[4]+'\t'+f[5]+'\t'+take1[0]+'\n')


    def n_map(self, yeast_id, domain): 

        """ maps different ids to nucleotides data; yeast_id=yeastID.txt domain = nucleotide.txt"""

        with open('id_nucleotide.txt', 'w') as  id_domain:
            with open(yeast_id, 'r') as fl:
                for fe in fl:
                    f = fe.split()
                    with open(domain,'r') as dp:
                        for d in dp:
                            d=d.split()
                            if len(f)>2:
                                if f[0]==d[0]:
                                    take=d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]
                                    if take>0:
                                        id_domain=open('id_nucleotide.txt', 'a')
                                        id_domain.write(take+'\n')


    def nucleotide_map(self, file1, file2):

        """ maps mutations to DNA-protein binding motifs; file2 = id_nucleotide.txt file1=mutation """
        
        with open('nucleotide_map.txt', 'w') as mp:
            with open(file1,'r') as f:
                for line in f:
                    line1 = line.split()
                    with open(file2, 'r') as f2:
                        for line2 in f2:
                            line2 = line2.split()
                            if line1[0] == line2[2]:
                                if int(line1[1]) >= int(line2[4]) and int(line1[1]) <= int(line2[5]):
                                    take = line2[0]+'\t'+line1[0]+'\t'+line2[4]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[6]+'\t'+'UniProt'
                                    if take > 0:
                                        mp = open('nucleotide_map.txt', 'w')
                                        summary = open('summary.txt', 'a+')
                                        mp.write(take+'\n')
                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[4]+'\t'+'Nucleotide-Binding'+'\t'+'UniProt'+'\n')


    def BioGrid(self):

        """ downloads BioGrid ids in raw form for further processing including mapping and web browsing
        WARNING: requires powerful machines to work with as its expensive to open in machines with low momey
        """

        url='http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,database(BioGrid)'
        request = urllib2.Request(url)
        request.add_header('User-Agent', 'Python %s')
        response = urllib2.urlopen(request)
        page = response.read(200000000)
        file= open('uniprot_bioGrid.txt','w')
        file.write(page)
        file.close()
    

    def preWeb(self, file, mutation ): 

        """ maps mutations to BioGrid ids; file ='uniprot_bioGrid.txt' ,  mutation =mutated_proteins.txt""" 

        with open('biog.txt', 'w') as out:
            with open(file, 'r') as fl:
                for f in fl:
                    f = f.rstrip().split()
                    if len(f) > 1:
                        i = f[1].split(';')
                        take = f[0]+'\t'+i[0]
                        take = take.split()
                        with open(mutation, 'r') as pro:
                            for p in pro:
                                p = p.split()
                                if take[0] == p[0]:
                                    take2 = take[0]+'\t'+take[1]+'\t'+'UniProt'
                                    if take2 > 0:
                                        out = open('biog.txt', 'a')
                                        out.write(take2+'\n')

    def Web(self, file1): 

        """ opens the BioGrid db in browser with as many tabs as mutated proteins, if possible; # file = 'biog.txt'""" 

        url = 'http://thebiogrid.org/'
        fl = open(file1, 'r')
        for f in fl:
            f = f.split()
            webbrowser.open(url + f[1])


    def pdb_c(self, file_1):

        """ Structure data filteration from UniProt; # file = uniprot_mod_raw.txt"""

        with open('pdb.txt', 'w') as stru:
            with open(file_1, 'r') as raw:
                for r in raw:
                    if not r.startswith('##'):
                        line = r.split()
                        if line[2] == 'Beta' and len(line[9]) > 1:
                            take = line[9].split('|')
                            take3 = line[0]+'\t'+line[2]+'\t'+line[4]+'\t'+line[5]+'\t'+take[1]
                            if len(take3)>0:
                                stru = open('pdb.txt', 'a')
                                stru.write(take3+'\n')
                                continue
                        if len(line) > 7 and line[2] == 'Helix' or line[2]=='Turn':
                            if len(line[8])>1:
                                tak = line[8].split('|')
                                tak3 = line[0]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+take[1]
                                if len(tak3)>0:
                                    stru = open('pdb.txt', 'a')
                                    stru.write(tak3+'\n')

    
    def mu_map(self):

        """ mutations from mutated proteins file mapped to the yeastID file"""

        with open('mutation_id.txt', 'w') as f:
            with open('mutation.txt') as mutation_file:
                for a in mutation_file:
                    a = a.split()
                    with open('yeastID.txt') as id:
                        for i in id:
                            i = i.split()
                            if len(i) > 2:
                                if a[0] == i[2]:
                                    take = i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+a[1]
                                    if take>0:
                                        f = open('mutation_id.txt', 'a')
                                        f.write(take+'\n')


    def pdb(self, file_pdb):

        """ This code maps mutations to the proteins structural regions;  file_pdb = pdb.txt"""

        with open('stru_mutation.txt', 'w') as s:
            with open(file_pdb, 'r') as raw:
                for i in raw:
                    i = i.split()
                    with open('mutation_id.txt') as mu: 
                        for m in mu:
                            m = m.split()
                            if i[0] == m[0]:
                                if int(i[2]) <= int(m[3]) and int(i[3]) >= int(m[3]):
                                    take = m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+i[1]+'\t'+i[2]+'\t'+m[3]+'\t'+i[3]+'\t'+i[4]+'\t'+'UniProt'
                                    if take > 0:
                                        s = open('stru_mutation.txt', 'a')
                                        summary = open('summary.txt', 'a+')
                                        s.write(take+'\n')
                                        summary.write(m[0]+'\t'+m[2]+'\t'+m[3]+'\t'+i[4]+'\t'+i[1]+'\t'+'UniProt'+'\n')


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#/////////////////// Annotated PTMs data from other resources than UniProt (know to play role in PPI and cross-talk) /////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#To get the mutational effects on PPI, PTM based crosstalk, Protein domains, we need to run the following data files from one local dict.; the data
#retrieved from PTMcode 2.0 and PTMfunc, for this reason. To run your own lis tagainst this program, all you need to do is to change the file name in
#the test varible and then you get to go, the output contains the pvalue of the GO terms effected by the mutations and also step wise protein output data
#to interpret your experiment."""
#This frame work contains the advance stage of mapping, where same one code can be used for the mapping to the different 
#PTM types, present at interface  and/or ppi.


def interface(file1, mutation):

    """PTM present at the interface of two proteins and known to play role in interation (Beltrao et al. Cell 2012)"""
    
    with open('interface_mutation.txt', 'w') as out:
        with open(file1, 'rU') as f:
            for l in f:
                line = l.split()
                if len(line) > 5:
                    take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[5]
                    take = take.split()
                    with open(mutation) as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[1] and m[1] == take[2]:
                                take2 = take[0]+'\t'+take[1]+'\t'+take[2]+'\t'+take[3]+'\t'+'PTMfunc'
                                fi = take2.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'Interface'+'\t'+'PTMfunc'+'\n')
                                            if take2 > 0:
                                                out = open('interface_mutation.txt', 'a')
                                                out.write(take2+'\n')
                                                
         

def ppi(file1,mutation):

    """ PTM present at the interface of two proteins and known to play role in interation (PTMfunc; Beltrao et al. Cell 2012)"""

    with open('ppi_mutation.txt', 'w') as out:
        with open(file1, 'rU') as f:
            for ls in f:
                line = ls.split()
                with open (mutation) as mu:
                    for m in mu:
                        m = m.split()
                        if len(line) > 7:
                            if m[0] == line[1] and m[1] == line[3]:
                                take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+'PTMfunc'
                                fi = take.split()
                                with open('yeastID.txt') as i:
                                    for di in i:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+'\t'+'PPI'+'\t'+'PTMfunc'+'\n')
                                            if take > 0:
                                                out = open('ppi_mutation.txt' , 'a+')
                                                out.write(take+'\n')
                                                continue
                            if m[0] == line[6] and m[1] == line[3]:
                                take2 = line[6]+'\t'+line[2]+'\t'+line[3]+'\t'+'PTMfunc'
                                fi = take2.split()
                                with open('yeastID.txt') as i:
                                    for di in i:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+'\t'+'PPI'+'\t'+'PTMfunc'+'\n')
                                            if take2 > 0:
                                                out = open('ppi_mutation.txt' , 'a+')
                                                out.write(take2+'\n')
                    
    
def withinPro(file2, mutation):

    """ PTMs (predicted to be) involved in the crosstalk within a given protein at baker's years (Minguez el 2012)"""

    with open('within_protein.txt', 'w') as file1:
        with open(file2, 'rU') as f:
            for l in f:
                line = l.split()
                if len(line)>19:
                    take = line[15]+'\t'+line[16]+'\t'+line[3]+'\t'+line[17]+'\t'+line[7]+'\t'+line[19]
                    take = take.split()
                    with open(mutation, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[1] and m[1]==take[3]:
                                take2 = take[0]+'\t'+take[1]+'\t'+take[2]+'\t'+take[3]+'\t'+'PTMcode'
                                fi=take2.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[3]+'\t'+fi[2]+'\t'+'WithinProtein'+'\t'+'PTMcode'+'\n')
                                            if take2 > 0:
                                                file1 = open('within_protein.txt', 'a')
                                                file1.write(take2+'\n')
                                                continue
                            if m[0] == take[1] and m[1] == take[5]:
                                take3 = take[0]+'\t'+take[1]+'\t'+take[4]+'\t'+take[5]+'\t'+'PTMcode'
                                fi = take3.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[3]+'\t'+fi[2]+'\t'+'WithinProtein'+'\t'+'PTMcode'+'\n')
                                            if take3 > 0:
                                                file1 = open('within_protein.txt', 'a')
                                                file1.write(take3+'\n')
                         

def betweenPro(fileb, mutation):

    """ PTMs (predicted to be) involved in the crosstalk in different proteins at baker's years (PTMcode 2.0; Minguez el 2012), 
    returns proteins present in both lists, only; file = ('sc_btw_proteins.txt', 'rU') mutation = open ('mutation') """

    with open('ptm_between_proteins.txt', 'w') as file1:
        with open(fileb, 'rU') as f:
            for l in f:
                line = l.split()
                if len(line)>20:
                    take = line[16]+'\t'+line[18]+'\t'+line[15]+'\t'+line[17]+'\t'+line[19]+'\t'+line[21]+'\t'+line[4]+'\t'+line[8]
                    take = take.split()
                    with open(mutation, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if m[0] == take[0] and m[1]==take[4]:
                                take2 = take[0]+'\t'+take[2]+'\t'+take[4]+'\t'+take[6]+'\t'+'PTMcode'
                                fi=take2.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take2 > 0:
                                                file1 = open('ptm_between_proteins.txt', 'a+')
                                                file1.write(take2+'\n')
                                                continue
                            if m[0] == take[1] and m[1] == take[5]:
                                take3 = take[1]+'\t'+take[3]+'\t'+take[5]+'\t'+take[7]+'\t'+'PTMcode'
                                fi=take3.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di=di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take3 > 0:
                                                file1 = open('ptm_between_proteins.txt', 'a+')
                                                file1.write(take3+'\n')

    
def hotspot(fileh, mutation):

    """ PTMs containing motifs in a close proximity are named hotspots (Beltrao et al. Cell 2012) 
    fileh = open('schotspot_updated.txt', 'rU') #mutation=open('mutation') """

    with open('hotspot.txt', 'w') as hotspot:
        with open(fileh, 'rU') as f:
            for l in f:
                line = l.split()
                with open(mutation, 'r') as mu:
                    for m in mu:
                        m = m.split()
                        if len(line) > 6:
                            if m[0] == line[2] and m[1] == line[3]:
                                take = line[1]+'\t'+line[2]+'\t'+line[3]+'\t'+line[5]+'\t'+line[6]+'\t'+'PTMfunc'
                                fi = take.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[1]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'HotSpot'+'\t'+'PTMFunc'+'\n')
                                            if take > 0:
                                                hotspot = open('hotspot.txt', 'a')
                                                hotspot.write(take+'\n')
                          

def sum_file_map():  

    """ reports all the results in a 'final-report' file """

    with open('final_report.txt', 'w') as x:
        with open('summary.txt') as fil1:
            for f in OrderedDict.fromkeys(fil1):
                if f > 0:
                    x.write(f)


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                               USEAGE       (Optional) 
#------------------------------------------------------------------------------------------------------------
#This usage strategy is optional, and a user can use above written codes in any convenient way as
#required by experiemental settings and data interpretation
#////////////////////////////////////////////////////////////////////////////////////////////////////////////



c = YGtPM()
wd = os.getcwd()
co=codon_from_protein_file()


def data(): 

    """ this function will download and clean required data to run ymap methods smoothly """

    start_time = time.time()
    try:
        dat = c.PTMdata()
    except IOError:
        pass
    try:
        cl = c.clean('uniprot_mod_raw.txt')
    except IOError:
        pass
    try:
        i = c.ID()
    except IOError:
        pass
    try:
        m = c.Map('yeastID.txt', 'PTMs.txt')
    except IOError:
        pass
    try:
        d = c.dclean('uniprot_mod_raw.txt')
    except IOError:
        pass
    try:
        dm = c.d_map('yeastID.txt', 'domains.txt')
    except IOError:
        pass
    try:
        ab = c.ab('uniprot_mod_raw.txt')
    except IOError:
            pass
    try:
        ii = c.id('bact.txt', 'yeast_id.txt')
    except IOError:
            pass
    try:
        bio=c.BioGrid()
    except IOError:
            pass
    try:
        c.pdb_c('uniprot_mod_raw.txt')
    except IOError:
        pass
    try:
        co.gff()
    except IOError:
        pass
    try:
        co.frmt('gff.txt')
    except IOError:
        pass
    try:
        co.id_map('yeastID.txt', 'frmt.txt')
    except IOError:
        pass
    try:
        c.nucleotide()
    except IOError:
        pass
    try:
        c.n_map('yeastID.txt', 'nucleotide.txt')
    except IOError:
        pass 
    try:
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_aceksites_interfaceRes_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_phosphosites_interfaceRes_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_ubisites_interfaceRessc_sc.txt', wd) 
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_acet_interactions.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'sc_btw_proteins.txt', wd) 
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_psites_interactions_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_ubi_interactions_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'sc_within_proteins.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'schotspot_updated.txt', wd)    
    except IOError:
        pass
    return "All required data downloaded in %s seconds" % (time.time() - start_time)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#//////////////////////////////// Following two codes are used for return the mutations at proteins level \\\\\\\\\\\\\\\\\\
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def mutation_types_file(): 

    """ mutation type and amino acid calculation where ref. and mutant base known """

    start_time = time.time()
    if not os.path.exists('mutated_proteins.txt'):
        raise StopIteration('because of missing mutation file with name "mutated_proteins.txt"')
    else:
        try:
            mutation_file('mutated_proteins.txt', 'd_id_map.txt')
        except IOError:
            pass
    return "Mutations with mutations types are available to map on functional entities"


def code(): 

    """ codon calculations where ref and mutant base unknown """

    start_time = time.time()
    try:
        co.cc_mapping('mutation1.txt', 'd_id_map.txt')
    except IOError:
            pass
    return "codons are available to map on functional entities"


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#////////////////////////////////// Following series of codes will return three files - mapped mutations, pvalue and biog.txt - for each type of data types \\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def ptm():

    """ PTMs mapping to mutations """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            a = c.ptm_map('mutation.txt', 'PTM_id_file.txt')
        except IOError:
            pass
        try:    
            p = c.enrich('mutated_proteins.txt')
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'mutated_proteins.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'PTMs')
            shutil.move(wd+"/"+'mutated_proteins.txt', wd+"/"+'PTMs')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs')
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'PTMs')
        except IOError:
            pass
    return "PTMs mapped in %s seconds" % (time.time() - start_time)


def domain():

    """ protein domain mapping """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            dom = c.dmap('mutation.txt', 'id_domain.txt')
        except IOError:
            pass
        try:
           p = c.enrich('domains_mapped.txt')  
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'domains_mapped.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'Domains')
        except IOError:
            pass
        try:
            shutil.move(wd+"/"+'domains_mapped.txt', wd+"/"+'Domains')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Domains')
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'Domains')
        except IOError:
            pass
    return "Domains mapped in %s seconds" % (time.time() - start_time)
    

def nucleo():

    """ DNA-protein binding motif mapping """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            c.nucleotide_map('mutation.txt', 'id_nucleotide.txt')
        except IOError:
            pass
        try:
           p = c.enrich('nucleotide_map.txt')  
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'nucleotide_map.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'Nucleotide_binding')
        except IOError:
            pass
        try:
            shutil.move(wd+"/"+'nucleotide_map.txt', wd+"/"+'Nucleotide_binding')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Nucleotide_binding')
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'Nucleotide_binding')
        except IOError:
            pass
    return "Nucleotide_binding domains mapped in %s seconds" % (time.time() - start_time)


def ab():

    """ active and binding site mapping """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            mm = c.mmap('sites_id.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('ab_mutation_file.txt')
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'ab_mutation_file.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'A-B-sites')
            shutil.move(wd+"/"+'ab_mutation_file.txt', wd+"/"+'A-B-sites')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'A-B-sites')
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'A-B-sites')
        except IOError:
            pass
    return "Active-Binding proteins sites mapped in %s seconds" % (time.time() - start_time)


def struc_map():

    """ structural regions mapping """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            c.mu_map()
        except IOError:
            pass
        try:
            pd = c.pdb('pdb.txt')
        except IOError:
            pass
        try:
            p = c.enrich('stru_mutation.txt')
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'stru_mutation.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'PDB')
            shutil.move(wd+"/"+'stru_mutation.txt', wd+"/"+'PDB')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PDB')
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'PDB')
        except IOError:
            pass
        return "Mutations are mapped to structural features in %s seconds" % (time.time() - start_time)

def intf():

    """ east = effective data which shows PTMs present at interface, ppi and 
        domain (hotspot) this analaysis could lead to an effective way to interpret
        user's mutational data on Yeast proteins from PTMfunc (also 3DID db) and PTMcode 2.0"""
    
    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            interface('3DID_aceksites_interfaceRes_sc.txt' ,'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('interface_mutation.txt')
        except IOError:
            pass
            #c.preWeb('uniprot_bioGrid.txt', 'interface_mutation.txt') 
        try:
            os.system("mkdir "+wd+"/"+'Interface')
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'acetylation')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'acetylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'acetylation')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'Domains')
        except IOError:
            pass
        try:
            interface('3DID_phosphosites_interfaceRes_sc.txt' ,'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('interface_mutation.txt')
        except IOError:
            pass
            #c.preWeb('uniprot_bioGrid.txt', 'interface_mutation.txt') 
        try:
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'Phosphorylation')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'Interface'+"/"+'Phosphorylation')
        except IOError:
            pass
        try:   
            interface('3DID_ubisites_interfaceRessc_sc.txt' ,'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('interface_mutation.txt')
        except IOError:
            pass
        #c.preWeb('uniprot_bioGrid.txt', 'interface_mutation.txt')
        try:
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'ubiquitination')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'Interface'+"/"+'ubiquitination')
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def pi():
    start_time = time.time()
    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            ppi('SC_acet_interactions.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            c.enrich('ppi_mutation.txt')
        except IOError:
            pass
        #c.preWeb('uniprot_bioGrid.txt', 'ppi_mutation.txt') 
        try:
            os.system("mkdir "+wd+"/"+'PPI')
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'acetylation')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'acetylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'acetylation')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PPI'+"/"+'acetylation')
        except IOError:
            pass
        try:
            ppi('SC_psites_interactions_sc.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('ppi_mutation.txt')
            #c.preWeb('uniprot_bioGrid.txt', 'ppi_mutation.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'Phosphorylation')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PPI'+"/"+'Phosphorylation')
        except IOError:
            pass
    
        try:
            ppi('SC_ubi_interactions_sc.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            c.enrich('ppi_mutation.txt')
        except IOError:
            pass
            #c.preWeb('uniprot_bioGrid.txt', 'ppi_mutation.txt') 
        try:
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'ubiquitination')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PPI'+"/"+'ubiquitination')
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def withP():
    start_time = time.time()
    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            withinPro('sc_within_proteins.txt','mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('within_protein.txt')
        except IOError:
            pass
        #c.preWeb('uniprot_bioGrid.txt', 'within_protein.txt')
        try:
            os.system("mkdir "+wd+"/"+'PTMs_within_Proteins')
            shutil.move(wd+"/"+'within_protein.txt', wd+"/"+'PTMs_within_Proteins')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_within_Proteins')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PTMs_within_Proteins')
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)

def betweenP():
    start_time = time.time()
    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            betweenPro('sc_btw_proteins.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('ptm_between_proteins.txt')
        except IOError:
            pass
        #c.preWeb('uniprot_bioGrid.txt', 'ptm_between_proteins.txt')
        try:
            os.system("mkdir "+wd+"/"+'PTMs_between_Proteins')
            shutil.move(wd+"/"+'ptm_between_proteins.txt', wd+"/"+'PTMs_between_Proteins')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_between_Proteins')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PTMs_between_Proteins')
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)


def hotS():
    start_time = time.time()
    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            hotspot('schotspot_updated.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('hotspot.txt')
        except IOError:
            pass
        #c.preWeb('uniprot_bioGrid.txt', 'hotspot.txt')
        try:
            os.system("mkdir "+wd+"/"+'PTMs_hotSpots')
            shutil.move(wd+"/"+'hotspot.txt', wd+"/"+'PTMs_hotSpots')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_hotSpots')
            #shutil.move(wd+"/"+'biog.txt', wd+"/"+'PTMs_hotSpots')
        except IOError:
            pass
        return "run time is %s seconds" % (time.time() - start_time)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Following two codes with perform all the codes on all the data /////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def uniprot_data(): 

    """ to perform all functions on UniProt(like ptm, domain and ab () functions) all together """ 

    try:
        ptm()
    except IOError:
        pass
    try:
        domain()
    except IOError:
        pass
    try:
        ab()
    except IOError:
        pass
    try:
        struc_map()
    except IOError:
        pass
    try:
        nucleo()
    except IOError:
        pass
    return "The Uniprot data is resolved into functional for interpretation"


def functional_data():

    """ to perform all functions on UniProt(like ptm, domain and ab () functions) all together """

    if not os.path.exists('mutation.txt'):
            raise StopIteration('because of missing mutation file')
    else:
        try:
            intf()
        except IOError:
            pass
        try:
            pi()
        except IOError:
            pass
        try:
            withP()
        except IOError:
            pass
        try:
            betweenP()
        except IOError:
            pass
        try:
            hotS()
        except IOError:
            pass
        return "The data from PTMcode and PTMfunc on PTMs functional biasedness is resolved into functional for interpretation"


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////  Final module of ymap package for executing all the modules \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


def ymap():

    """ returns all the results of all the codes pf yMap """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        try:
            uniprot_data()
        except IOError:
            pass
        try:
            functional_data()
        except IOError:
            pass
        try:
            sum_file_map()
        except IOError:
            pass
        try:
            y = (time.time() - start_time)
            os.makedirs('yMap-results'+str(y))
        except IOError:
            pass
        try:    
            p = c.enrich('final_report.txt')
        except IOError:
            pass
        try:
            c.preWeb('uniprot_bioGrid.txt', 'final_report.txt')
        except IOError:
            pass
        try:
            shutil.move(wd+"/"+'PTMs', wd+"/"+'yMap-results'+str(y))
            shutil.move('Domains', wd+"/"+'yMap-results'+str(y))
            shutil.move('Nucleotide_binding',wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'A-B-sites', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'PDB', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'Interface', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'PPI', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'PTMs_within_Proteins', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'PTMs_between_Proteins',wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'PTMs_hotSpots',wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'mutation.txt', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'final_report.txt', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'biog.txt', wd+"/"+'yMap-results'+str(y))
            shutil.move(wd+"/"+'summary.txt', wd+"/"+'yMap-results'+str(y))
        except IOError: 
            pass
        return "All functional data is ready in about %s seconds" % (time.time() - start_time)


def web(): 
    """ NOTE: to use the following function change to dir to respective folder to run web based analysis """   
    if not os.path.exists('biog.txt'):
        raise StopIteration('because of missing mutation file')
    else:
        c.Web('biog.txt')
    return "Web is ready for networks exploration of mutated proteins"




