SIGMAS = []
AVGS = []
for file in files:
    A = n.load(file)
    SIGMAS.append(A['SIGMA'])
    AVGS.append(A['AVG'])
    freqs = A['freqs']*1e3
figure(31,figsize=(6,5))
clf()
ch = n.argwhere(n.logical_and(freqs>120,freqs<180))
for i,AVG in enumerate(AVGS):
    plot(freqs[ch],n.real(AVG[ch])/n.max(n.real(AVGS)),label='day %d'%i)
SIGMAS=n.array(SIGMAS) 
SIGMA = n.sqrt(n.mean(SIGMAS**2,axis=0))[ch].squeeze()/n.max(n.real(AVGS))
AVG = n.mean(AVGS,axis=0)[ch].squeeze()/n.max(n.real(AVGS))
F = freqs[ch].squeeze()
y1 = AVG-SIGMA/2
y2 = AVG+SIGMA/2
fill_between(F,y1=n.real(y1),y2=n.real(y2),
    facecolor='0.3',
    edgecolor='0.3',
    zorder=20,
    alpha=0.9)
plot(freqs[ch],n.real(n.mean(AVGS,axis=0)[ch]),'k')
ylabel('$\mathcal{R}e(V_{8,12})$  [counts]')
xlabel('frequency [MHz]')
legend(loc='lower right',ncol=4,columnspacing=0.5,numpoints=4,prop={'size':10})
subplots_adjust(left=0.15,bottom=0.15)
