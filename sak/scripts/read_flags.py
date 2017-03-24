import numpy as np, sys, matplotlib.pyplot as plt

files = sys.argv[1:]

c = np.zeros((112),dtype='int')
grid = np.zeros((len(files),112),dtype='int')
for i,f in enumerate(files):
    with open(f) as f: lines = f.readlines()
    jd,bas = lines[0].rstrip().split('    ')
    jd,bas = int(jd), map(int,bas.split(','))
    print jd,bas
    for ba in bas:
        c[ba]+=1 #XXX must be an array manipulation way to do this
        grid[i,ba] = 1

plt.imshow(grid,aspect='auto',interpolation='nearest')
plt.xlabel('Antenna number')
plt.ylabel('JD')
plt.suptitle('red=bad, blue=good')
plt.show()
plt.close()

plt.step(range(112),c,where='mid')
plt.xlim(-0.5,112)
plt.xlabel('Antenna number')
plt.ylabel('Badness count over Season')
plt.show()
        
