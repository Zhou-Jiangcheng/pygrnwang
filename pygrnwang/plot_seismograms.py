import matplotlib.pyplot as plt
import numpy as np


def plot_seismograms(seismograms, srate, ylabel="u (m)", cut_length=None):
    if cut_length is None:
        length = len(seismograms[0])
    else:
        length = cut_length * srate
    t = np.arange(0, round(length / srate), 1 / srate)
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(t, seismograms[0][:length], color="red", linewidth=0.2)
    plt.subplot(3, 1, 2)
    plt.plot(t, seismograms[1][:length], color="green", linewidth=0.2)
    plt.subplot(3, 1, 3)
    plt.plot(t, seismograms[2][:length], color="blue", linewidth=0.2)
    plt.xlabel("t (s)")
    plt.ylabel(ylabel)
    plt.show()
