alnfile = '/scratch/dmoi/datasets/covid_data/msa_0501/msa_0501.fasta'
count = 0
with open(alnfile , 'r') as infile:
    with open( alnfile+'.cleaned.fasta' , 'w') as outfile:
        for l in infile:
            if '>' in l:
                lclean= l.replace(' ', '_').replace('/' , '|')
                outfile.write( '>'+ str(count) + lclean[1:] )
                count+=1
            else:
                outfile.write(l)
#write to philip informative
from Bio import SeqIO
count +=1
    with open( alnfile+'.cleaned.fasta' , 'w') as outfile:
        
def write_phy(fasta , count):
    total_len =
with open( alnfile+'.cleaned.fasta', "r") as input_handle:
    with open(alnfile + 'phy' , "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "fasta")
        count = SeqIO.write(sequences, output_handle, "phylip")
print("Converted %i records" % count)
