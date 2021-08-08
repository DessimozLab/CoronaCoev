alnfile = '/scratch/dmoi/datasets/covid_data/30_may/mmsa_2021-06-01/2021-06-01_masked.fa'


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