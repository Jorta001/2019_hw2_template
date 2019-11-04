#!/usr/bin/env python

import os, gzip, itertools, sys, re, collections, pprint

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhim$
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

stdoutOrigin=sys.stdout
sys.stdout = open("Se_genes.txt", "w")

# QUESTION 1

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)

    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        print(seqname, " first 10 bases are ", seqstring[0:10])

sys.stdout.close()
sys.stdout=stdoutOrigin

fhand=open('Se_genes.txt')
count=0
for line in fhand:
        count=count+1
print ("Question 1a: The total number of genes in Salmonella enterica is", count)

stdoutOrigin=sys.stdout
sys.stdout = open("Mt_genes.txt", "w")


for l in gzip.open(file2,"rt").readlines():
        if not l.startswith(">") and not l.startswith("...skipping one line"): print(l)

sys.stdout.close()
sys.stdout=stdoutOrigin

with open('Mt_length.txt') as infile:
    lines=0
    words=0
    characters=0
    for line in infile:
        wordslist=line.split()
        lines=lines+1
        words=words+len(wordslist)
        characters += sum(len(word) for word in wordslist)
print("Question 2b:The total length of gene sequences for Mycobacterium tuberculosis is", characters)

# Question 3

from string import ascii_lowercase
from collections import Counter

with open('Se_length.txt') as f:
    print("Question 3a: The frequencies of nitrogenous bases in Salmonella Enterica")
    print Counter(letter for line in f 
                  for letter in line.lower() 
                  if letter in ascii_lowercase)

with open('Mt_length.txt') as f:
    print("Question 3b The frequencies of nitrogenous bases in Mycobacterium tubercolusis")
    print Counter(letter for line in f

...skipping one line
                  if letter in ascii_lowercase)

# QUESTION 4

# this reponse is mostly based off of this stackoverflow response https://stackoverflow.com/questions/46385576/python-how-to-count-frequency-of-triple-nucleotide

stdoutOrigin=sys.stdout
sys.stdout = open("Se_Codons.txt", "w")

with open('Se_length.txt', "r") as dna_sequence:

	def parse_sequence(dna_sequence):
		codons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG', 'TAT', 'TAC', 'TAA', 'TAG', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'GCT','GCC', 'GCA', 'GCG', 'TGT', 'TGC', 'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']
		if len(dna_sequence) % 3 == 0:
                	for i in range(0,len(dna_sequence),3):
                        	codons.append((dna_sequence[i:i + 3]))
		return codons

	def frequency(dna_sequence):
		parsed = parse_sequence(dna_sequence)
		return ct.Counter(parsed)

sys.stdout.close()
sys.stdout=stdoutOrigin
