import matplotlib.pyplot as plt
import numpy as np


def plot_seismograms(seismograms, srate, length_in_s=None):
    if length_in_s is None:
        length = len(seismograms[0][:])
    else:
        length = length_in_s * srate
    t = np.arange(0, round(length / srate), 1 / srate)
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(t, seismograms[0][:length], color="red", linewidth=0.2)
    plt.subplot(3, 1, 2)
    plt.plot(t, seismograms[1][:length], color="green", linewidth=0.2)
    plt.subplot(3, 1, 3)
    plt.plot(t, seismograms[2][:length], color="blue", linewidth=0.2)
    plt.xlabel("t [s]")
    plt.ylabel("u [m]")
    plt.show()
