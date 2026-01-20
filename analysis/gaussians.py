import matplotlib.pyplot as plt
import numpy as np


def pF(t, tmax, sigma, C, offset=0.):
    return C * np.exp(-1/2 * (t-tmax)**2 / sigma**2) + offset
   

def cumulInteg(t, f):
    N = len(t)
    intervals = N-1
    h = t[1] - t[0]
    out = np.zeros(shape=(N,))
    extra = (t[1:] - t[:-1]) / 2 + t[:-1]
    intval = np.zeros(shape=(intervals+N,))
    intval[0::2] = t
    intval[1::2] = extra 
    for i in range(0, int((len(intval)-1)/2), 2):
        out[i] = f[i] + 4*f[i+1] + f[i+2]
    return np.cumsum(out) * h/3

def AnalyticalComparison():
    N = 199 # must be uneven for Simpson
    t = np.linspace(0, 30, N)

    C0 = 35
    f = pF(t, 15, 4, 3 , 0) + pF(t, 5, 3, 4, 0)
    I = cumulInteg(t, f) + C0
    samples = [5, 10, 15, 20, 25, 30]
    Iint = np.interp(samples, t, I)
    return Iint

if __name__ == '__main__':

    N = 999 # must be uneven for Simpson
    t = np.linspace(0, 30, N)

    C0 = 35
    f = pF(t, 15, 4, 3 , 0) + pF(t, 5, 3, 4, 0)
    I = cumulInteg(t, f) + C0
    samples = [5, 10, 15, 20, 25, 30]
    Iint = np.interp(samples, t, I)
    print("\n".join([f"t :  {samples[i]:2d}, I : {Iint[i]}" for i in range(len(samples))]))


    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 6))
    fig.suptitle("Prozessfunktion")
    ax1.plot(t, f, c="b", label="Quellterm")
    ax1.set_title("Quellfunktion")
    ax1.set_xlabel("t [s]") 
    ax1.set_ylabel(r"$\partial_t T$ [°C/s]")
    ax2.plot(t, I, c="k", label="Integ. Quelle")
    ax2.set_title("Integrierte Temperatur")
    ax2.set_xlabel("t [s]")
    ax2.set_ylabel(" T [°C]")
    ax1.legend()
    plt.savefig("ana.png")
    #plt.show()
