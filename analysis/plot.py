import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gaussians import AnalyticalComparison


df = pd.read_excel("project.ods", engine="odf")

t = df["t[s]"].to_numpy()
top05 = df["top 05D"].to_numpy()
top1 = df["top1D"].to_numpy()
top5 = df["top5D"].to_numpy()

tops = [top05, top1, top5]

simp05 = df["simp05D"].to_numpy()
simp1 = df["simp1D"].to_numpy()
simp5 = df["simp5D"].to_numpy()

simps = [simp05, simp1, simp5]

l = [tops, simps]

c = ["b", "g", "o"]
label = [r"$D=0.5$", r"$D=1.$", r"$D=5.$"]

fig, axes = plt.subplots(
    ncols=2, figsize=(8, 6), sharex=True, sharey=True
)


for i in range(len(axes)):
    for j in range(len(c)):
        axes[i].plot(
            t, l[i][j],
            ls="--", marker="x", label=label[j]
        )
    # end loop
        
    axes[i].plot(
            t, AnalyticalComparison(),
            ls="--", marker="x", c="k", label="ana comp."
        )
    axes[i].set_xlabel(r"$t$ [s]")
    axes[i].set_ylabel(r"$T$ [°C]")
    axes[i].legend()

axes[0].set_title("CPU top")
axes[1].set_title("CPU simple")

fig.suptitle("Temperaturmaxima in dem Volumen")
plt.savefig("Tcomp.png")
#plt.show()
