import os,random
import numpy as np

output_dataset_file = 'dataset.txt'
output_query_file = 'query.txt'

dimension = 5
numPoints = 1000

if os.path.exists(output_dataset_file):
	os.remove(output_dataset_file)

dim_wise_factorial = []

for j in range(dimension):
	dim_wise_factorial.append(np.random.choice(numPoints,numPoints,replace = False)/numPoints)

# print(dim_wise_factorial)

with open(output_dataset_file,'w+') as outFile:
	s = str(dimension) +' '+ str(numPoints) +'\n'
	outFile.write(s)
	for i in range(numPoints):
		s = ''
		for j in range(dimension):
			s = s + str(round(dim_wise_factorial[j][i],3)) + ' '
		s = s + '\n'
		# print(s)
		outFile.write(s)


with open(output_query_file,'w+') as outFile:
	s = str(dimension) + '\n'
	for j in range(dimension):
		s = s + str(round(random.uniform(-1,1),1)) + ' '
	s = s+'\n'
	outFile.write(s)	
