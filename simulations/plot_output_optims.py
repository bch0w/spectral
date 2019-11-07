import os
import numpy as np
import matplotlib.pyplot as plt

fids = ["./cc/output.optim_cc", "./mtm/output.optim_mtm"]
for fid in fids:
    iterations, steplens, misfits = [], [], []
    i = 0
    iters = []
    lines = open(fid, "r").readlines()
    for line in lines[2:]:
        line = line.strip().split()
        # Each iteration will have an integer to represent the iter number
        if len(line) == 3:
            iteration = line[0]
            iterations.append(int(iteration))
            steplens.append(float(line[1]))
            misfits.append(float(line[2]))
            i+=1 
            iters.append(i)
        # Each trial step length will follow and will carry the same iteration
        elif len(line) == 2:
            iterations.append(int(iteration))
            steplens.append(float(line[0]))
            misfits.append(float(line[1]))
            i+=1 
            iters.append(i)
        elif len(line) == 0:
            continue
        else:
            print(line)
            print("invalid line length encountered in output.optim")
    
    plt.plot(iters, misfits, 'o-', label=os.path.basename(fid))
    plt.legend()
    plt.title('Optimization convergence MTM vs CC')
    plt.xlabel('Step lengths')
    plt.ylabel('Pyatoa Misfits')
    plt.grid(True)

plt.savefig('30event_mtmcc.png')
