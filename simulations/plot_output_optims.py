import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.6
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['axes.linewidth'] = 2


def plot_by_itr():
    """
    Plot by iteration, ignoring intermediate step count information
    """
    fids = ["./output.optim_cc", "./output.optim_mtm", "./output.optim_both"]
    for i, fid in enumerate(fids):
        itr_misfits = []
        lines = open(fid, "r").readlines()
        for line in lines[2:]:
            line = line.strip().split()
            # Each iteration will have an integer to represent the iter number
            if len(line) == 3:
                itr_misfits.append(float(line[2]))
        
        # Plot iteration and misfits
        itrs = np.linspace(1, len(itr_misfits), len(itr_misfits))
        plt.plot(itrs, itr_misfits, 'o-', label=os.path.basename(fid), 
                 zorder=10)
        plt.legend()
        plt.title('Optimization convergence; 30event 79rec')
        plt.xlabel('Iteration')
        plt.ylabel('Pyatoa Misfits')
        plt.grid(True)

    plt.savefig('30event_mtmcc_itr.png')
    plt.close("all")

def plot_by_step():
    """
    Plot by step count, plot iterations as vertical annotated lines
    """
    fids = ["./output.optim_cc", "./output.optim_mtm", "./output.optim_both"]
    for i, fid in enumerate(fids):
        step, itr = 0, 0
        steps, misfits = [], []
        lines = open(fid, "r").readlines()
        for line in lines[2:]:
            line = line.strip().split()
            # Each iteration will have an integer to represent the iter number
            if len(line) == 3:
                misfits.append(float(line[2]))
                itr += 1 
                step += 1
                # Mark where each iteration starts
                if i == 0:
                    plt.axvline(step, linestyle='--', c='k', zorder=2)
                    plt.annotate(s="m{:0>2}".format(itr-1), xy=(step, 0.85),
                                 fontsize=8)
                steps.append(step)
            # Each trial step length will follow and will carry the same iteration
            elif len(line) == 2:
                misfits.append(float(line[1]))
                step += 1 
                steps.append(step)
            elif len(line) == 0:
                continue
            else:
                print(line)
                print("invalid line length encountered in output.optim")
        
        plt.plot(steps, misfits, 'o-', label=os.path.basename(fid), zorder=10)
        plt.legend()
        plt.title('Optimization convergence; 30events 79rec')
        plt.xlabel('Step lengths')
        plt.ylabel('Pyatoa Misfits')
        plt.grid(True)

    plt.savefig('30event_mtmcc_step.png')
    plt.close("all")


if __name__ == "__main__":
    plot_by_step()
    plot_by_itr()
