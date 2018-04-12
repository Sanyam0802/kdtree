import os,random
import numpy as np
import sys, math

output_dataset_file = 'dataset.txt'
output_queries = 'query.txt'
output_query_file = 'query_final.txt'

if len(sys.argv) != 4:
	print("Invalid Format")
	print("python file_generator.py <dimension> <numPoints> <numQueryPoints>")
	exit()

dimension = int(sys.argv[1])
numPoints = int(sys.argv[2])
numQueryPoints = int(sys.argv[3])

if os.path.exists(output_dataset_file):
	os.remove(output_dataset_file)

if os.path.exists(output_queries):
	os.remove(output_queries)

if os.path.exists(output_query_file):
	os.remove(output_query_file)

# dim_wise_factorial = []

# for j in range(dimension):
# 	dim_wise_factorial.append(np.random.choice(numPoints,numPoints,replace = False)/numPoints)
# print(len(dim_wise_factorial[0]))

# print(dim_wise_factorial)

with open(output_dataset_file,'w+') as outFile:
	s = str(dimension) +' '+ str(numPoints) +'\n'
	outFile.write(s)
	round_range = math.ceil(math.log10(numPoints))
	# print(round_range)
	for i in range(numPoints):
		s = ''
		for j in range(dimension):
			# s = s + str(round(dim_wise_factorial[j][i],round_range)) + ' '
			s = s + str(round(random.uniform(0,1),round_range)) + ' '
		s = s + '\n'
		# print(s)
		outFile.write(s)

with open(output_queries,'w+') as outFile:
	s = str(dimension) +' '+ str(numQueryPoints) +'\n'
	outFile.write(s)
	for i in range(numQueryPoints):
		s = ''
		for j in range(dimension):
			s = s + str(round(random.random(),3)) + ' '
		s = s + '\n'
		# print(s)
		outFile.write(s)


with open(output_query_file,'w+') as outFile:
	s = str(dimension) + '\n'
	for j in range(dimension):
		s = s + str(round(random.uniform(-1,1),1)) + ' '
	s = s+'\n'
	outFile.write(s)	
