import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

coverage = pd.read_csv('./pair_end_exp/taxa_extract.coverage', sep='\t', names=['genome', 'loci', 'depth'])

loci = coverage['loci'].values
depth = coverage['depth'].values

graph_loci = []
graph_depth = []

# Where going to bin loci locations and depths into groups of 5000. 
# From these bins of 5000 where going to calculate the mean of bin depths
# and then take the start loci of the bin.

# start_pos = 0
# stop_pos = 5000
# iterations = len(loci) / stop_pos
# print iterations


# for n in range(iterations):

#     if n == iterations - 1:
#         print "fetching the last bits of depth and loci"
#         print n
#     else:
#         print "running"


plt.scatter(loci, depth)

plt.title('Genome Coverage Plot: Treponema')

plt.xlabel('Loci')
plt.ylabel('Depth of Coverage')

plt.savefig('coverage_test.png')