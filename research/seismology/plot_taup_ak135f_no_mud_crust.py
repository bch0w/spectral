import matplotlib.pyplot as plt

with open(f"taup_models/ak135f_no_mud.nd", "r") as f:
    lines = f.readlines()[:7]

depths = []
vss = []
vpp = []
for line in lines:
    try:
        depth, vp, vs, rho, qk, qm = line.strip().split()
    except ValueError:
        continue
    
    depths.append(-1*float(depth))
    vss.append(float(vs))
    vpp.append(float(vp))


plt.plot(vss, depths, "ro-", label="Vs")
plt.plot(vpp, depths, "bo-", label="Vp")
plt.legend()
plt.xlabel("Velocity [km/s]")
plt.ylabel("Depth [km]")
plt.title("TauP ak135f_no_mud")
plt.show()
    


