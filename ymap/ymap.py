#!/usr/bin/env python
#
# year of release ! 2016
#
# With the exponential growth of Posttranslational modifications (PTMs) 
# data and and lack of characterisation of all the PTM-types. Its important
# to undersand properly the functions and experimental relevence of PTMs by 
# creating the tools that facilitate the PTMs based analyses. 
# And to understand the importance of PTMs in Yeast genome, its important 
# to make it easier to map experimental mutations to PTM positional data. 
# it's also important and relevent to translate genetic abrretions to understand 
# the phenotype. 
# We architect a python (yMap) library to help users to understand which parts of
# mutated proteins are affected during the yeast experimentation.
# This facilitation not only would help bioligists to interpret their data 
# efficiently but also gives freedom to save time by mapping data to mutations 
# easily 
#
# The yMap program is a python based fast and robust automated method to map 
# large yeast variants to proteins post-translational modifications, proteins domains,
# proteins-DNA binding domains, proteins structural regions, proteins active and 
# binding sites, proteins networks visualisation. 
# For Usage see README file
#
# Dependencies:
# Orange Bioinformatics 
# see README file 
#

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals
try:
    from builtins import next
    from builtins import str
    from builtins import range
    from builtins import object
    from builtins import bytes
except ImportError:
    pass
import os
import math
import zipfile
from itertools import groupby
import shutil
import time
import urllib
from collections import OrderedDict
import webbrowser
try:
    import Orange
except ImportError:
    import Orange3
from orangecontrib.bio import go    
from six.moves import range
ontology = go.Ontology()
annotations = go.Annotations("sgd", ontology=ontology)
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen



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


def mutation_file(mutation, d_id):
        """ defines the mutation types; either Non-Synonmous or Stop Codon"""
        with open('mutation.txt', 'wb') as t:   
            with open(mutation, 'rU') as mut: 
                for m in mut:
                    m = m.rstrip().split()
                    with open(d_id,'rU') as id:    
                        for i in id:
                            i = i.rstrip().split()
                            if not m[0].startswith('c'.upper()):
                                if len(m) != 5  or not m[0].startswith('c'.lower()):
                                    raise StopIteration('Please correct the format of input mutation file')
                                else:
                                    if m[4] == i[2]:
                                        take = m[4]+'\t'+m[0]+'\t'+i[3]+'\t'+m[1]+'\t'+i[4]+'\t'+m[2]+'\t'+m[3]+'\t'+i[5]
                                        take1= take.rstrip().split()       
                                        with open('gff.txt', 'rU') as orf:     
                                            linee = orf.readlines()[23078:]
                                            up = (x[1] for x in groupby(linee, lambda line: line[0] == ">")) 
                                            for head in up:
                                                head = next(head)[1:].strip()
                                                seq = "".join(s.strip() for s in next(up))
                                                if head == take1[1] and take1[0] == i[2] and take1[7] == '-':   
                                                    cod = 1 + (int(take1[4])-int(take1[3]))                       
                                                    cc = math.ceil(int(cod)/float(3))
                                                    c = str(cc).split('.')                         
                                                    cn = int(c[0])-1    
                                                    sli_n = seq[int(take1[2]):int(take1[4])]                  
                                                    rev_sli_n = revcomp(sli_n, reverse=True, complement=True)  
                                                    sli_m_n = sli_n[:int(-cod)]+take1[6]+sli_n[int(-cod)+1:] 
                                                    rev_sli_m_n = revcomp(sli_m_n, reverse=True, complement=True)   
                                                    wild_type_rev_n = translate_dna(rev_sli_n)                
                                                    mut_type_n = translate_dna(rev_sli_m_n)
                                                    try:
                                                        if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] == '_':
                                                            pic = take1[0]+'\t'+str(c[0])+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Stop' +'\t'+take1[1]+'\t'+take1[3]
                                                            if pic > str(0): 
                                                                t = open('mutation.txt', 'a')
                                                                t.write(pic+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                    try:
                                                        if wild_type_rev_n[cn] != mut_type_n[cn] and mut_type_n[cn] != '_':
                                                            pic = take1[0]+'\t'+str(c[0])+'\t'+wild_type_rev_n[cn]+'\t'+mut_type_n[cn]+'\t'+'Non-Synonymous' +'\t'+take1[1]+'\t'+take1[3]                                                                                                                    
                                                            if pic > str(0):
                                                                t = open('mutation.txt', 'a+')
                                                                t.write(pic+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                if head == take1[1] and take1[0]==i[2] and take1[7] == '+':
                                                    code = int(take1[3])-int(take1[2])
                                                    code1 = 1 + (int(take1[3])-int(take1[2])) 
                                                    cce = math.ceil(int(code1)/float(3)) 
                                                    ce = str(cce).split('.') 
                                                    cp = int(ce[0])-1                  
                                                    pos = int(take1[2]) - 1                               
                                                    sli_p = seq[int(pos):int(take1[4])]                   
                                                    sli_m_p = sli_p[:int(code)]+take1[6]+sli_p[int(code)+1:]  
                                                    wild_type_p = translate_dna(sli_p)                    
                                                    mut_type_p = translate_dna(sli_m_p)
                                                    try:                   
                                                        if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp] != '_': 
                                                            pick = take1[0]+'\t'+str(ce[0])+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Non-Synonymous'+'\t'+take1[1]+'\t'+take1[3]
                                                            if pick > str(0):
                                                                with open('mutation.txt', 'a+') as t:
                                                                    t.write(pick+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic1+'\n')
                                                        continue
                                                    try:
                                                        if wild_type_p[cp] != mut_type_p[cp] and mut_type_p[cp]=='_':
                                                            pick = take1[0]+'\t'+str(ce[0])+'\t'+wild_type_p[cp]+'\t'+mut_type_p[cp]+'\t'+'Stop' +'\t'+take1[1]+'\t'+take1[3]
                                                            if pick > str(0):
                                                                with open('mutation.txt', 'a+') as t:
                                                                    t.write(pick+'\n')
                                                    except IndexError as e:
                                                        pic1 =  take1[0]+ '\t'+ 'Error:'+'\t'+ str(e)
                                                        t = open('mutation.txt', 'a+')
                                                        t.write(pic1+'\n')


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# //////////////////                     UniProt data                 /////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class YGtPM(object):

    
    def gff(self):
        """ The genomic coordinates downloaded in gff formate for further processing to calculate mutated codons, if not
              available, see next method"""
        rsponse = urlopen("http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff")
        page = rsponse.read()
        file = open('gff.txt','wb')  
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
                            with open('frmt.txt','a') as file4:
                                file4.write(result2+'\n')


    def id_map(self, file_id, frmt):        
        with open('d_id_map.txt', 'w') as file2:
            with open(file_id, 'r') as file_id_name:
                for line in file_id_name:
                    line=line.split()
                    with open(frmt, 'r') as fr:
                        for f in fr:
                            f=f.split()
                            if len(line)>2:
                                if line[1]==f[0]:
                                    result= line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+f[1]+'\t'+f[2]+'\t'+f[3]
                                    if result > str(0):
                                        with open('d_id_map.txt', 'a') as file2:
                                            file2.write(result+'\n')


    def pTMdata(self):

        """Downloads UpiProt data as a raw txt file (uniprot_mod_raw.txt)"""
        rsponse = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=gff&columns=id,feature(MODIFIED%20RESIDUE)')
        page = rsponse.read()
        fil = open('uniprot_mod_raw.txt','wb')
        fil.write(page)
        fil.close()
    
    
    def clean(self, UniProt_file):              
        
        """ cleans file uniprot_mod_raw.txt into a tab separated PTMs.txt
        """

        with open('PTMs.txt', 'w') as out:
            with open(UniProt_file,'rU') as UniProt_file_name:
                for l in UniProt_file_name:
                    if not l.startswith('##'):
                        line = l.split()
                        if line[2] == 'Lipidation':
                            lll = line[0]+'\t'+line[4]+'\t'+line[8]
                            ll = lll.split()
                            ll = ll[2].split('=')
                            p = line[0]+'\t'+line[4]+'\t'+ll[1]
                            if p > str(0):
                                out = open('PTMs.txt', 'a')
                                out.write(p+'\n')
                                continue
                        if line[2] == 'Glycosylation':
                            ggg = line[0]+'\t'+line[4]+'\t'+line[8]
                            gg = ggg.split()
                            gg =  gg[2].split('=')
                            p1 = line[0]+'\t'+line[4]+'\t'+gg[1]
                            if p1 > str(0):
                                out = open('PTMs.txt', 'a+')
                                out.write(p1+'\n')
                                continue
                        if line[2] == 'Modified':
                            mmm = line[0]+'\t'+line[4]+'\t'+line[9]
                            mm = mmm.split()
                            mm = mm[2].split('=')
                            mm = mm[1].split(';')
                            p2 = line[0]+'\t'+line[4]+'\t'+mm[0]
                            if p2 > str(0):
                                out = open('PTMs.txt', 'a+')
                                out.write(p2+'\n')
                                continue
                        if line[2] == 'Cross-link': #ubiquitination
                            ccc = line[0]+'\t'+line[4]+'\t'+line[8]
                            cc = ccc.split()
                            cc = cc[2].split('=')
                            p3 = line[0]+'\t'+line[4]+'\t'+cc[1]
                            if p3 > str(0):
                                with open('PTMs.txt', 'a+') as out:
                                    out.write(p3+'\n')


    def iD(self):

        """ This method retrieves the different ID types for maping """
        rsponse = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,genes(OLN),%2Cgenes(PREFERRED)')
        page1 = rsponse.read()
        file_1 =open('yeastID.txt','wb')  
        file_1.write(page1)
        file_1.close()
        
            
    def pmap(self, file_id, file_PTMs):          

        """ if proteins ids are not SDG or uniprot or common names, this method maps the ids 
        """
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
                                    if result3 > str(0):
                                        file3 = open('PTM_id_file.txt', 'a')
                                        file3.write(result3+'\n')
                                             

    def ptm_map(self, mutation_file, PTM_id_file):

        """ This method maps the overlap between mutated codons from previous method to the PTM sites"""
        summary = open('summary.txt', 'w')
        with open('mutated_proteins.txt', 'w') as file5:
            with open(mutation_file, 'rU') as mutation_file:
                for line in mutation_file:
                    line = line.split()
                    with open(PTM_id_file, 'r') as file_PTMs:
                        for line1 in file_PTMs:
                            line1 = line1.split()
                            if line[0] == line1[2] and line[1] == line1[3]:
                                take = line1[0]+'\t'+line1[1]+'\t'+line[0]+'\t'+line[1]+'\t'+line1[4]+'\t'+'UniProt'
                                if take > str(0):
                                    file5 = open('mutated_proteins.txt', 'a')
                                    summary = open('summary.txt', 'a')
                                    file5.write(take+'\n')
                                    summary.write(line1[0]+'\t'+line[0]+'\t'+line[1]+'\t'+line1[4]+'\t'+'PTMs'+'\t'+'UniProt'+'\n')


    def dclean(self, uniprot_mod_raw):  

        """domain data needed to be filters from UniProt file, before mapping domains"""

        with open('domains.txt', 'w') as domain:
            with open(uniprot_mod_raw, 'rU') as raw:
                for a in raw:
                    if not a.startswith('##'):
                        a = a.split('=')
                        a1 = a[0].split()
                        if a1[2] == 'Domain':
                            if len(a) == 2:
                                a2 = a[1].rstrip()
                                take = a1[0]+'\t'+a1[3]+'\t'+a1[4]+'\t'+a2+'\n'
                                if take > str(0):
                                    with open('domains.txt', 'a') as domain:
                                        domain.write(take)
                                        continue
                            if len(a) == 4:
                                a3 = a[1].rstrip().split(';')
                                a4 = a[3].rstrip().split('|')
                            if len(a4) > 1:
                                take2 = a1[0]+'\t'+a1[3]+'\t'+a1[4]+'\t'+a3[0]+'\t'+a4[1]+'\n'
                                if take2 > str(0):
                                    with open('domains.txt', 'a+') as domain:
                                        domain.write(take2)
                            if len(a4) == 1:
                                take3 = a1[0]+'\t'+a1[3]+'\t'+a1[4]+'\t'+a4[0]+'\n'
                                if take3 > str(0):
                                    with open('domains.txt', 'a+') as domain:
                                        domain.write(take3)
                
    
    def d_map(self, yeast_id, domain):

        """ maps the different proteins ids to domains"""
        
        with open('id_domain.txt','w') as id_domain:
            with open(yeast_id, 'rU') as fl:
                for f in fl:
                    f = f.split()
                    with open(domain,'r') as dp:
                        for d in dp:
                            d = d.split()
                            if len(f) > 2 and f[0] == d[0]:
                                if len(d) == 4:
                                    take = d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\n'
                                    if take > str(0):
                                        with open('id_domain.txt','a') as id_domain:
                                            id_domain.write(take)
                                if len(d) == 5:
                                    take1 = d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]+'\n'
                                    if take1 > str(0):
                                        with open('id_domain.txt','a+') as id_domain:
                                            id_domain.write(take1)
                                if len(d) == 6:
                                    take2 = d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]+'\t'+d[5]+'\n'
                                    if take2 > str(0):
                                        with open('id_domain.txt','a+') as id_domain:
                                            id_domain.write(take2)
                                if len(d) == 7:
                                    take3 = d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]+'\t'+d[5]+'\t'+d[6]+'\n'
                                    if take3 > str(0):
                                        with open('id_domain.txt','a+') as id_domain:
                                            id_domain.write(take3)
                                if len(d) == 8:
                                    take4 = d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]+'\t'+d[5]+'\t'+d[6]+'\t'+d[7]+'\n'
                                    if take4 > str(0):
                                        with open('id_domain.txt','a+') as id_domain:
                                            id_domain.write(take4)


    def dmap(self, file1, file2):

        """ maps mutations to the yeast domains"""
        
        with open('domains_mapped.txt', 'w') as mp:
            with open('summary.txt', 'a+') as summary:
                with open(file1,'rU') as f:
                    for line in f:
                        line1=line.split()
                        with open(file2, 'r') as f2:
                            for line2 in f2:
                                line2 = line2.split()
                                if line1[0] == line2[2]:
                                    try:
                                        if line1[1] == 'Error:':
                                            with open('domains_mapped.txt', 'a') as mp:
                                                mp.write("input file contains error position for" +line1[0]+ "this protein"+'\n')
                                                continue
                                        if int(line1[1]) >= int(line2[3]) and int(line1[1]) <= int(line2[4]):
                                            if len(line2) == 6:
                                                take = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+'UniProt'+'\n'
                                                if take > str(0):
                                                    with open('domains_mapped.txt', 'a+') as mp:
                                                        mp.write(take)
                                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+'domain'+'\t'+'UniProt'+'\n')
                                            if  len(line2) == 7:
                                                take1 = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+line2[6]+'UniProt'+'\n'
                                                if take1 > str(0):
                                                    with open('domains_mapped.txt', 'a+') as mp:
                                                        mp.write(take1)
                                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+line2[6]+'domain'+'\t'+'UniProt'+'\n')
                                            if  len(line2) == 8:
                                                take2 = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'UniProt'+'\n'
                                                if take2 > str(0):
                                                    with open('domains_mapped.txt', 'a+') as mp:
                                                        mp.write(take2)
                                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'domain'+'\t'+'UniProt'+'\n')
                                            if  len(line2) == 9:
                                                take3 = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'\t'+line2[8]+'UniProt'+'\n'
                                                if take3 > str(0):
                                                    with open('domains_mapped.txt', 'a+') as mp:
                                                        mp.write(take3)
                                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'\t'+line2[8]+'domain'+'\t'+'UniProt'+'\n')
                                            if  len(line2) == 10:
                                                take4 = line2[0]+'\t'+line1[0]+'\t'+line2[3]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'\t'+line2[8]+'\t'+line2[9]+'UniProt'+'\n'
                                                if take4 > str(0):
                                                    with open('domains_mapped.txt', 'a+') as mp:
                                                        mp.write(take4)
                                                        summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[5]+'\t'+line2[6]+'\t'+line2[7]+'\t'+line2[8]+'\t'+line2[9]+'domain'+'\t'+'UniProt'+'\n')
                                    except IndexError:
                                        pass
    
    def enrich(self, file1):

        """ This method performed enrichment analysis of mutated proteins and
        return the p value of functional enrichment of mutated proteins functional regions/residues; 
        see main text for how pvalue is calculated"""

        with open('pvalue.txt','w') as out:
            with open(file1, 'rU') as f:
                for line in f:
                    line=line.split()
                    res = annotations.get_enriched_terms(line)
                    if len(res) == 0:
                        out = open('pvalue.txt','a')
                        out.write('No enrichment found')
                    else:
                        for go_id, (genes, p_value, ref) in list(res.items()):
                            if p_value < 0.002:
                                take = ontology[go_id].name, p_value, ''.join(genes)
                                with open('pvalue.txt','a') as out:
                                    out.write(str(take))
                                    continue
                            if p_value < 0.005:
                                take1 = ontology[go_id].name, p_value, ''.join(genes)
                                with open('pvalue.txt','a+') as out:
                                    out.write(str(take1))
                                    continue
                            if p_value < 0.02:
                                take2 = ontology[go_id].name, p_value, ''.join(genes)
                                with open('pvalue.txt','a+') as out:
                                    out.write(str(take2))
                                    continue
                            if p_value < 0.05:
                                take3 = ontology[go_id].name, p_value, ''.join(genes)
                                with open('pvalue.txt','a+') as out:
                                    out.write(str(take3))


    def ab(self, file_raw): 

        """Prepares raw Uniprot data for yeast active and binding sites mutation analysis"""

        with open('bact.txt','w') as file2:
            with open(file_raw, 'rU') as d:
                for f in d:
                    if not f.startswith('##'):
                        f = f.split()
                        if f[2] == 'Active':
                            take = f[0]+'\t'+f[2]+'\t'+f[4]
                            if take > str(0):
                                with open('bact.txt','a') as file2:
                                    file2.write(take+'\n')
                        if f[2] == 'Binding':
                            take2 = f[0]+'\t'+f[2]+'\t'+f[4]
                            if take2 > str(0):
                                with open('bact.txt','a+') as file2:
                                    file2.write(take2+'\n')
                            
                             
    def id(self, act, yeast_id): 

        """ maps proteins ids to active and binding sites containing proteins"""

        with open('sites_id.txt', 'w') as file_id:
            with open(act, 'rU') as a:
                for a in a:
                    a = a.split()
                    with open('yeastID.txt')as id:
                        for i in id:
                            i = i.split()
                            if len(i) > 2:
                                if a[0] == i[0]:
                                    take = i[2]+'\t'+i[1]+'\t'+i[0]+'\t'+a[1]+'\t'+a[2]
                                    if take > str(0):
                                        with open('sites_id.txt', 'a') as file_id:
                                            file_id.write(take+'\n')


    def mmap(self, file_sites, mutation):

        """ maps mutations to proteins ab (active and binding sites) """ 

        with open('ab_mutation_file.txt', 'w') as out: 
            with open(file_sites, 'rU') as s:
                for a in s:
                    a = a.split()
                    with open(mutation, 'r') as mu:
                        for m in mu:
                            m = m.split()
                            if a[0] == m[0] and a[4] == m[1]:
                                take = a[2]+'\t'+ a[3]+'\t'+m[1]+'\t'+'UniProt'
                                if take > str(0):
                                    with open('ab_mutation_file.txt', 'a') as out: 
                                        summary = open('summary.txt', 'a+')
                                        out.write(take+'\n')
                                        summary.write(a[2]+'\t'+a[0]+'\t'+m[1]+'\t'+ a[3]+'\t'+'Active/Binding site'+'\t'+'UniProt'+'\n')


    def nucleotide(self):

        """ prepares the UniProt data for the nucleotide motifs mapping to mutations """

        with open('nucleotide.txt', 'w') as t:
            with open('uniprot_mod_raw.txt', 'rU') as file_raw:
                for fi in file_raw:
                    if not fi.startswith('##'):
                        f = fi.split()
                        if f[2] == 'Nucleotide' and len(f) > 8:
                            take = f[0]+'\t'+f[2]+'\t'+f[4]+'\t'+f[5]+'\t'+f[9]
                            take1 = take.split()
                            take1 = take1[4].split(';')
                            if take > str(0):
                                with open('nucleotide.txt', 'a') as t:
                                    t.write(f[0]+'\t'+f[2]+'\t'+f[4]+'\t'+f[5]+'\t'+take1[0]+'\n')


    def n_map(self, yeast_id, domain): 

        """ maps different proteins ids to nucleotides data """

        with open('id_nucleotide.txt', 'w') as  id_domain:
            with open(yeast_id, 'rU') as fl:
                for fe in fl:
                    f = fe.split()
                    with open(domain,'r') as dp:
                        for d in dp:
                            d=d.split()
                            if len(f)>2:
                                if f[0]==d[0]:
                                    take=d[0]+'\t'+f[1]+'\t'+f[2]+'\t'+d[1]+'\t'+d[2]+'\t'+d[3]+'\t'+d[4]
                                    if take > str(0):
                                        with open('id_nucleotide.txt', 'a') as  id_domain:
                                            id_domain.write(take+'\n')


    def nucleotide_map(self, file1, file2):

        """ maps mutations to protein-nucleotide binding motifs """
        
        with open('nucleotide_map.txt', 'w') as mp:
            with open(file1,'rU') as f:
                for line in f:
                    line1 = line.split()
                    with open(file2, 'r') as f2:
                        for line2 in f2:
                            line2 = line2.split()
                            if line1[0] == line2[2]:
                                try:
                                    if line1[1] == 'Error:':
                                        with open('nucleotide_map.txt', 'a') as mp:
                                            mp.write("input file contains error position for" + line1[0]+"protein"+'\n')
                                            continue
                                    if int(line1[1]) >= int(line2[4]) and int(line1[1]) <= int(line2[5]):
                                        take = line2[0]+'\t'+line1[0]+'\t'+line2[4]+'\t'+line1[1]+'\t'+line2[4]+'\t'+line2[6]+'\t'+'UniProt'
                                        if take > str(0):
                                            with open('nucleotide_map.txt', 'a') as mp:
                                                summary = open('summary.txt', 'a+')
                                                mp.write(take+'\n')
                                                summary.write(line2[0]+'\t'+line1[0]+'\t'+line1[1]+'\t'+line2[4]+'\t'+'Nucleotide-Binding'+'\t'+'UniProt'+'\n')
                                except IndexError:
                                    pass


    def bioGrid(self):

        """ Downloads BioGrid ids of yeast proteins from UniProt for further processing including mapping and web browsing
        WARNING: requires powerful machines to work with as its expensive to open in machines with low memory
        """
        response = urlopen('http://www.uniprot.org/uniprot/?query=yeast&fil=organism%3A%22Saccharomyces%20cerevisiae%20(strain%20ATCC%20204508%20%2F%20S288c)%20(Baker%27s%20yeast)%20%5B559292%5D%22&sort=score&format=tab&columns=id,database(BioGrid)')
        page = response.read()
        file1 = open('uniprot_bioGrid.txt','wb')
        file1.write(page)
        file1.close()
    

    def preWeb(self, file1, mutation ): 

        """ maps mutations to BioGrid ids """ 

        with open('biog.txt', 'w') as out:
            with open(file1, 'rU') as fl:
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
                                    if take2 > str(0):
                                        with open('biog.txt', 'a') as out:
                                            out.write(take2+'\n')

    def bweb(self, file1): 

        """ opens the BioGrid db in browser with as many tabs as mutated proteins""" 

        url = 'http://thebiogrid.org/'
        fl = open(file1, 'rU')
        for f in fl:
            f = f.split()
            webbrowser.open(url + f[1])


    def pdb_c(self, file_1):

        """ Structure data filtration from UniProt"""

        with open('pdb.txt', 'w') as stru:
            with open(file_1, 'rU') as raw:
                for r in raw:
                    if not r.startswith('##'):
                        line = r.split()
                        if line[2] == 'Beta' and len(line[9]) > 1:
                            take = line[9].split('|')
                            take3 = line[0]+'\t'+line[2]+'\t'+line[4]+'\t'+line[5]+'\t'+take[1]
                            if take3 > str(0):
                                with open('pdb.txt', 'a') as stru:
                                    stru.write(take3+'\n')
                                    continue
                        if len(line) > 7 and line[2] == 'Helix' or line[2]=='Turn':
                            if len(line[8])>1:
                                tak = line[8].split('|')
                                tak3 = line[0]+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+take[1]
                                if tak3 > str(0):
                                    with open('pdb.txt', 'a+') as stru:
                                        stru.write(tak3+'\n')

    
    def mu_map(self):

        """ mutations proteins mapped to the yeastID file"""

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
                                    if take > str(0):
                                        with open('mutation_id.txt', 'a') as f:
                                            f.write(take+'\n')


    def pdb(self, file_pdb):

        """ This code maps mutations to the proteins structural regions"""

        with open('stru_mutation.txt', 'w') as s:
            with open(file_pdb, 'rU') as raw:
                for i in raw:
                    i = i.split()
                    with open('mutation_id.txt') as mu: 
                        for m in mu:
                            m = m.split()
                            if i[0] == m[0]:
                                try:
                                    if m[3] == 'Error:':
                                        with open('stru_mutation.txt', 'a') as s:
                                            s.write("input file contains error position for" +m[2]+ "protein"+'\n')
                                            continue
                                    if int(i[2]) <= int(m[3]) and int(i[3]) >= int(m[3]):
                                        take = m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+i[1]+'\t'+i[2]+'\t'+m[3]+'\t'+i[3]+'\t'+i[4]+'\t'+'UniProt'
                                        if take > str(0):
                                            with open('stru_mutation.txt', 'a+') as s:
                                                summary = open('summary.txt', 'a+')
                                                s.write(take+'\n')
                                                summary.write(m[0]+'\t'+m[2]+'\t'+m[3]+'\t'+i[4]+'\t'+i[1]+'\t'+'UniProt'+'\n')
                                except IndexError:
                                    pass


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

    """PTM present at the interface of two proteins and known to play role in interaction (Beltrao et al. Cell 2012)"""
    
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
                                            if take2 > str(0):
                                                with open('interface_mutation.txt', 'a') as out:
                                                    out.write(take2+'\n')
                                                
         

def ppi(file1,mutation):

    """ PTM present at the interface of two proteins and known to play role in interaction (PTMfunc; Beltrao et al. Cell 2012)"""

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
                                            if take > str(0):
                                                with open('ppi_mutation.txt', 'a') as out:
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
                                            if take2 > str(0):
                                                with open('ppi_mutation.txt', 'a+') as out:
                                                    out.write(take2+'\n')
                    
    
def withinPro(file2, mutation):

    """ PTMs (predicted) involved in the crosstalk within a given protein at baker's years (Minguez el 2012)"""

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
                                            if take2 > str(0):
                                                with open('within_protein.txt', 'a') as file1:
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
                                            if take3 > str(0):
                                                with open('within_protein.txt', 'a+') as file1:
                                                    file1.write(take3+'\n')
                         

def betweenPro(fileb, mutation):

    """ PTMs (predicted) involved in the crosstalk in different proteins at baker's years (PTMcode 2.0; Minguez el 2012) """

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
                                fi = take2.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take2 > str(0):
                                                with open('ptm_between_proteins.txt', 'a') as file1:
                                                    file1.write(take2+'\n')
                                                    continue
                            if m[0] == take[1] and m[1] == take[5]:
                                take3 = take[1]+'\t'+take[3]+'\t'+take[5]+'\t'+take[7]+'\t'+'PTMcode'
                                fi=take3.split()
                                with open('yeastID.txt') as id:
                                    for di in id:
                                        di = di.split()
                                        if len(di) > 2 and di[2] == fi[0]:
                                            summary = open('summary.txt', 'a+')
                                            summary.write(di[0]+'\t'+di[2]+'\t'+fi[2]+'\t'+fi[3]+'\t'+'BetweenProteins'+'\t'+'PTMcode'+'\n')
                                            if take3 > str(0):
                                                with open('ptm_between_proteins.txt', 'a+') as file1:
                                                    file1.write(take3+'\n')

    
def hotspot(fileh, mutation):

    """ PTMs containing motifs in a close proximity are named hotspots (Beltrao et al. Cell 2012)"""

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
                                            if take > str(0):
                                                with open('hotspot.txt', 'a') as hotspot:
                                                    hotspot.write(take+'\n')
                          

def sum_file_map():  

    """ reports all the results in a 'final-report' file """

    with open('final_report.txt', 'w') as x:
        with open('summary.txt') as fil1:
            for fi in OrderedDict.fromkeys(fil1):
                x.write(fi)


#////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                               USEAGE       (Optional) 
#------------------------------------------------------------------------------------------------------------
#This usage strategy is optional, and a user can use above written codes in any convenient way as
#required by experiemental settings and data interpretation (see README for proper use)
#////////////////////////////////////////////////////////////////////////////////////////////////////////////


c = YGtPM()
wd = os.getcwd()

def data(): 

    """ this function will download and clean required data to run ymap methods smoothly """

    start_time = time.time()
    try:
        dat = c.pTMdata()
    except IOError:
        pass
    try:
        cl = c.clean('uniprot_mod_raw.txt')
    except IOError:
        pass
    try:
        i = c.iD()
    except IOError:
        pass
    try:
        m = c.pmap('yeastID.txt', 'PTMs.txt')
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
        bio=c.bioGrid()
    except IOError:
            pass
    try:
        c.pdb_c('uniprot_mod_raw.txt')
    except IOError:
        pass
    try:
        c.gff()
    except IOError:
        pass
    try:
        c.frmt('gff.txt')
    except IOError:
        pass
    try:
        c.id_map('yeastID.txt', 'frmt.txt')
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
        z = zipfile.ZipFile(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'sc_btw_proteins.txt.zip', 'r')
        z.extractall()
    except IOError:
        pass 
    try:
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_aceksites_interfaceRes_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_phosphosites_interfaceRes_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'3DID_ubisites_interfaceRessc_sc.txt', wd) 
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_acet_interactions.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_psites_interactions_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'SC_ubi_interactions_sc.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'sc_within_proteins.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'schotspot_updated.txt', wd)
        shutil.copy2(wd+'/'+'data'+'/'+'PTMcode+PTMfunc_data'+'/'+'sc_btw_proteins.txt', wd)    
    except IOError:
        pass
    return "All required data downloaded in %s seconds" % (time.time() - start_time)


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#//////////////////////////////// Following two codes are used for return the mutations at proteins level \\\\\\\\\\\\\\\\\\
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def mutation_types_file(): 

    """ mutation type and amino acid change calculation where ref. and mutant base known """

    start_time = time.time()
    try:
        mutation_file("mutated_proteins.txt", 'd_id_map.txt')
    except IOError:
        pass
    return "Mutations with mutations types are available to map on functional entities"


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
        try:
            os.system("mkdir "+wd+"/"+'Interface')
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'acetylation')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'acetylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'acetylation')
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
        try:
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'Phosphorylation')
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
        try:
            os.system("mkdir "+wd+"/"+'Interface'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'interface_mutation.txt', wd+"/"+'Interface'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'Interface'+"/"+'ubiquitination')
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
        try:
            os.system("mkdir "+wd+"/"+'PPI')
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'acetylation')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'acetylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'acetylation')
        except IOError:
            pass
        try:
            ppi('SC_psites_interactions_sc.txt', 'mutation.txt')
        except IOError:
            pass
        try:
            p = c.enrich('ppi_mutation.txt')
        except IOError:
            pass
        try:
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'Phosphorylation')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'Phosphorylation')
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
        try:
            os.system("mkdir "+wd+"/"+'PPI'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'ppi_mutation.txt', wd+"/"+'PPI'+"/"+'ubiquitination')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PPI'+"/"+'ubiquitination')
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
        try:
            os.system("mkdir "+wd+"/"+'PTMs_within_Proteins')
            shutil.move(wd+"/"+'within_protein.txt', wd+"/"+'PTMs_within_Proteins')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_within_Proteins')
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
        try:
            os.system("mkdir "+wd+"/"+'PTMs_between_Proteins')
            shutil.move(wd+"/"+'ptm_between_proteins.txt', wd+"/"+'PTMs_between_Proteins')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_between_Proteins')
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
        try:
            os.system("mkdir "+wd+"/"+'PTMs_hotSpots')
            shutil.move(wd+"/"+'hotspot.txt', wd+"/"+'PTMs_hotSpots')
            shutil.move(wd+"/"+'pvalue.txt', wd+"/"+'PTMs_hotSpots')
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


def ymap_genes():

    """ returns all the results of all the codes of yMap; starting from genetics coordinates of proteins """

    start_time = time.time()
    if not os.path.exists('mutation.txt'):
        try:
            mutation_types_file()
        except IOError:
            pass
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
            os.remove(wd+"/"+'mutation_id.txt')
        except IOError: 
            pass
        return "All functional data from genomic coordinates is ready in about %s seconds" % (time.time() - start_time)


def ymap_proteins():

    """ returns all the results of all the codes of yMap; starting from proteins level mutation positions """

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
            os.remove(wd+"/"+'mutation_id.txt')
        except IOError: 
            pass
        return "All functional data from proteins mutation-positions is ready in about %s seconds" % (time.time() - start_time)


def web(): 
    """ NOTE: to use the following function change to dir to respective folder to run web based analysis """   
    if not os.path.exists('biog.txt'):
        os.chdir(input('specify biog.txt path:'))
        c.bweb('biog.txt')
    return "Web is ready for networks exploration of mutated proteins"

def path():
    "Path to the BioGrid ids path for visualisation"
    try:
        os.chdir(raw_input("paste here path to biog.txt file:"))
    except IOError:
        pass
    return "you need to provide path/to/biog.txt"

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-d', '--data', help='downloads required data to run yMap successfully')
    parser.add_argument('-g','--ymap_genes', help='performs the yMap on genes level mutation positions')
    parser.add_argument('-p', '--ymap_proteins', help='performs the yMap on proteins level mutation positions')
    parser.add_argument('-w', '--web', help='generates BioGrid web pages for interactome visualisation; paste the path to biog.txt file')

    args = parser.parse_args()
    if args.data:
        try:
            data()
        except IOError:
            pass
    elif args.ymap_genes:
        try:
            ymap_genes()
        except IOError:
            pass
    elif args.ymap_proteins:
        try:
            ymap_proteins()
        except IOError:
            pass
    elif args.web: 
        os.chdir(input('specify biog.txt path:'))   # specify biog.txt path:'/yMap-results78.50792193412781'
        try:
            web()
        except IOError:
            pass
    else:
        print ("to run a function seek help")

