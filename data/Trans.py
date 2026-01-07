import os, sys

# filename = sys.argv[1]
filename = "RS_21_12_e.matrix"

with open(filename, 'r') as din:
    content = din.read().split('\n')
    x ,y = content[0].split(' ')
    m = content[1:int(x)+1]
    print(y+" "+x)
    for i in range(int(y)):
        for j in range(int(x)):
            print(m[j][i], end='')
        print("")