import sys
files = sys.argv[1:]
file_tuples = []
for file in files: file_tuples.append(file.split('.'))
S = sorted(file_tuples, key=lambda time: time[2])
for tup in S: print '.'.join(tup),
