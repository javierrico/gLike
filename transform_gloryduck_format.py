import numpy as np

lines = []
with open('gLike_outputfile.txt') as inputfile:
    lines = inputfile.readlines()
    lines = np.array(lines)

data = []
for line in np.arange(lines.shape[0]):
    data.append(lines[line].split(' '))
data = np.array(data)
output_file = open("gloryduck_format.txt","w")
for i in np.arange(data.shape[1]-1):
     for j in np.arange(data.shape[0]):
         # write to file in gloryduck format
         output_file.write("{} ".format(data[j][i]))
     output_file.write("\n")
output_file.close()

