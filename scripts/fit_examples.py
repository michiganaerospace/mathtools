import numpy as np
from mathtools.fit import *
from vanity import *


if __name__ == '__main__':

    # Get our plots on!
    setup_plotting()

    # Create some noisy data.
    t = np.linspace(0,5*np.pi, 300)
    y = np.sin(2*np.pi/5*t) + 0.2 * np.random.randn(len(t))
    f = Fit(t,y, 15, reg_coefs=[0,1e-3,1e-4])

    figure(100)
    plot(t, y, 'o', markerfacecolor=pomegranate, markeredgecolor=pomegranate,\
         alpha=0.6)
    grid(True)
    xlabel('Time (seconds)')
    ylabel('Amplitude (volts)')
    savefig('docs/images/noisy_sine.png')

    
    figure(200)
    plot(t, y, 'o', markerfacecolor=pomegranate, markeredgecolor=pomegranate,\
         alpha=0.6)
    hold(True)
    plot(f.results.x, f.results.y, color=belize_hole, linewidth=2)
    grid(True)
    xlabel('Time (seconds)')
    ylabel('Amplitude (volts)')
    savefig('docs/images/noisy_sine_fit.png')




