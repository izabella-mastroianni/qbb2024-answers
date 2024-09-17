# week 1 

# step 1.1 How many 100bp reads are needed to sequence a 1Mbp genome to 3x coverage?
1Mbp*3 = 3Mbp
3Mbp/100bp = 30,000 reads

# 1.4
49988 bp in the genome have 0x coverage
It generally fits the Poisson expectations, but not very closely.
The normal distribution fits the data slightly better. 

# 1.5
116 bp in the genome have 0x coverage
The Poisson generally fits the data, and is closer than the 3x.
The normal distribution is still a better fit than Poisson, and also closer than 3x.

# 1.6
14 bp in the genome have 0x coverage
The Poisson expectations and normal distribution are good fits for the data and closely follow the distribution. 

# exercise 2

# 2.4
# run in unix
`dot -Tpng 2_1.dot > ex2_digraph.png`

# 2.5 
ATTGATTCTTATTCATTT

# 2.6
To more accurately reconstruct the sequence of the genome, you could do long read sequencing. Having longer reads makes it quicker and more accurate by making it easier to reconstruct repetetive parts of the genome. 

