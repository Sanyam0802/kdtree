import os,random

output_dataset_file = 'dataset.txt'
output_query_file = 'query.txt'

dimension = 3
numPoints = 10

if os.path.exists(output_dataset_file):
	os.remove(output_dataset_file)

with open(output_dataset_file,'w+') as outFile:
	s = str(dimension) +' '+ str(numPoints) +'\n'
	outFile.write(s)
	for i in range(numPoints):
		s = ''
		for j in range(dimension):
			s = s + str(round(random.uniform(-1,1),1)) + ' '
		s = s + '\n'
		outFile.write(s)


with open(output_query_file,'w+') as outFile:
	s = str(dimension) + '\n'
	for j in range(dimension):
		s = s + str(round(random.uniform(-1,1),1)) + ' '
	s = s+'\n'
	outFile.write(s)	