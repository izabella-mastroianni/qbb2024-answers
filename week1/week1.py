#!/usr/bin/env python3

import numpy
import scipy



# exercise 1, step 1.2

# for 3x coverage
genome_size = 1000000 #1Mbp genome size
read_size = 100 #100bp
coverage = 3 #3x coverage 
# 1Mbp * 3x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_size)


# array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads): 
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1


# save coverages as txt for R plotting
numpy.savetxt("genome_coverage.3x.txt", genome_coverage)


# range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage+1))

# get poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# get normal pdf at each of these
normal_estimates = scipy.stats.norm.pdf(xs,  
    numpy.mean(genome_coverage),
    numpy.std(genome_coverage))


# 1.5, 10x 
# for 10x coverage
genome_size = 1000000 #1Mbp genome size
read_size = 100 #100bp
coverage = 10 #10x coverage 
# 1Mbp * 10x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_size)


# array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads): 
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1


# save coverages as txt for R plotting
numpy.savetxt("genome_coverage.10x.txt", genome_coverage)


# range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage+1))

# get poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# get normal pdf at each of these
normal_estimates = scipy.stats.norm.pdf(xs,  
    numpy.mean(genome_coverage),
    numpy.std(genome_coverage))

# 1.5, 30x 
# for 30x coverage
genome_size = 1000000 #1Mbp genome size
read_size = 100 #100bp
coverage = 30 #30x coverage 
# 1Mbp * 30x coverage / 100bp reads
num_reads = int((genome_size * coverage)/read_size)


# array to keep track of the coverage at each position in the genome
genome_coverage = numpy.zeros(genome_size, int)

for _ in range(num_reads): 
    start_pos = numpy.random.randint(0, genome_size - read_size + 1)
    end_pos = start_pos + read_size
    genome_coverage[start_pos:end_pos] += 1


# save coverages as txt for R plotting
numpy.savetxt("genome_coverage.30x.txt", genome_coverage)


# range of coverages observed
max_coverage = max(genome_coverage)
xs = list(range(0, max_coverage+1))

# get poisson pmf at each of these
poisson_estimates = scipy.stats.poisson.pmf(xs, coverage)

# get normal pdf at each of these
normal_estimates = scipy.stats.norm.pdf(xs,  
    numpy.mean(genome_coverage),
    numpy.std(genome_coverage))


# exercise 2
# 2.1
reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']
graph = {}
k = 3

for read in reads: 
    for i in range(len(read) -k):
        kmer1 = read[i: i+k]
        kmer2 = read[i+1: i+1+k]
        graph.setdefault((kmer1,kmer2), 0)
        graph[(kmer1, kmer2)] += 1
print("digraph{")
for left, right in graph:
    print(f"{left} -> {right}")
print("}")

