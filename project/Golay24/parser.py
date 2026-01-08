import os,sys
import numpy

data_table = []
allr = []
cmp_idx = -int(sys.argv[1], 10)
for f in os.listdir('.'):
    if f.startswith("final_log_"):
        with open(f, 'r') as din:
            c = din.read().split('\n')[:-1]
        r = [float(s.split(', ')[cmp_idx]) for s in c[::2]]
        allr.extend([x for x in zip(r, c[::2], c[1::2])])
        r.sort()
        data_table.append((f, r))

with open("data_result.csv", 'w') as dout:
    for f, l in data_table:
        dout.write("%s,%s\n"%(f,",".join(["%lf"%s for s in l])))

allr.sort(key = lambda x: x[0])
min_r = allr[0][0]
print("min: %lf in %d case, var %lf"%(min_r, len(allr), numpy.var(numpy.array([x[0] for x in allr]))))
for a,b,c in allr:
    if a==min_r:
        print(b)
        print(c)