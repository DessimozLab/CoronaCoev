

alnfile = '/scratch/dmoi/datasets/covid_data/msa_0501/msa_0501.fasta'
idset = set()
count = 0
with open(alnfile , 'r') as infile:
    with open( alnfile+'.cleaned.fasta' , 'w') as outfile:
        for l in infile:
            if '>' in l:
                lclean= l.replace(' ', '_').replace('/' , '|')
                if lclean in idset:
                    outfile.write( lclean +str(count))
                    count+=1
                else:
                    outfile.write( lclean )
                    idset.add(lclean)
            else:
                outfile.write(l)
