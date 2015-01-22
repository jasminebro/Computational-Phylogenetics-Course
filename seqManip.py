"""
Program name:seqManip.py
Author:Jasmine Brown
python 2.7
"""
#define a DNA sequence-sequence from assignment page on github
sequence="aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggctg"
#This determines the length of my sequence 
print "The length of my sequence is", len(sequence)


#Now I need to create and store the RNA equivalent of the sequence, then print to screen.
#try to use .replace method to change out a letter for the rna of the string
RNAseq=sequence.replace("t","u") 
print RNAseq

#Create and store the reverse complement of your sequence, then print to the screen.
#I replace the original nucleotide with a capital one and replace that with its complement
Complement= sequence.replace("a","A").replace("A","T").replace("t","a").replace("g","G").replace("G","C").replace("c","g")
ReverseComplement=Complement[::-1] 
#now I print everything to the screen in lower case 
print "This is my reverse complement",Reversecomplement.lower()

# Extract the bases corresponding to the 13rd and 14th codons from the
#sequence, then print them to the screen. I used the index numbers for the codons and made a variable equal tot hat index
Codon13= sequence[36:39]
print Codon13
Codon14= sequence[39:42]
print Codon14

#Create a function to translate the nucleotide sequence to amino acids using the vertebrate mitochondrial genetic code 
###AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#Starts = ----------------------------------MM----------------------------
#Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

#Make a dictionary that holds the values for the amino acid along with its corresponding codon
AminoAcids={'FTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
            'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGC':'W','TGG':'W',
            'CTT':'T','CTG':'T','CTA':'T','CTG':'T','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCA':'A',
            'GAT':'D','GAC':'D','GAC':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

#this sequence has 2 nucleotides trimmed off			
newsequence="aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"
#I made all the nucleotides uppercase so it can match everything properly in the dictionary.
trimmedsequence=newsequence.upper()

#this determines my new sequence length
seqlength=len(trimmedsequence) 
print seqlength
#this created a blank variable to store the protein sequence as the loop goes through the codon list 
protein = ""

#the for loop starts at 0 and goes to the end of the trimmed sequence in increments of 3.    
for start in range(0,618,3):
#here we define a variable for the codons. Each codon comes from the trimmed sequence and starts at the zero index and counts the codons by 3 and the start position moves evey time it loops
    codon = trimmedsequence[start:start+3]
#a variable for the amino acid names is set. the variable goes to the dictionary and finds the set of 3 nucleotides from the codon variable and gets the AA that corresponds
    aa = AminoAcids.get(codon)
#I'm not sure if I typed some of the codons in correctly, but I would get a NONE value for some so if I got a NONE value I told it to keep going and translate the sequence
    if aa==None:
        continue
#this is adding the amino acids to the empty protein file above the loop as the loop goes through the trimmed sequence 
    protein = protein + aa
#this line prints my final protein 
print protein
