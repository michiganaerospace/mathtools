import numpy as np
from mathtools.fit import *
from mathtools.vanity import *


if __name__ == '__main__':

    # Get our plots on!
    setup_plotting()

    # Create some noisy data.
    t = np.linspace(0,5*np.pi, 300)
    t_new = np.linspace(-pi, 6*np.pi, 400)
    noise_amplitude = 0.1
    y = np.sin(2*np.pi/5*t) + noise_amplitude * np.random.randn(len(t))
    reg_coefs = [0, 0, 0]
    f = Fit(t, 15, reg_coefs=reg_coefs)

    # Perform the fit to the data.
    fit = f.fit(y)

    # Resample the fit to a new domain (including new points!)
    rfit = f.resample(t_new)


    # Make some plots to illustrate how our fitting works.
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
    plot(fit.x, fit.y, '.-', color=belize_hole, linewidth=2)
    grid(True)
    xlabel('Time (seconds)')
    ylabel('Amplitude (volts)')
    savefig('docs/images/noisy_sine_fit.png')


    figure(300)
    plot(t, y, 'o', markerfacecolor=pomegranate, markeredgecolor=pomegranate,\
         alpha=0.6)
    hold(True)
    plot(rfit.x, rfit.y, '.-', color=belize_hole, linewidth=2)
    grid(True)
    xlabel('Time (seconds)')
    ylabel('Amplitude (volts)')
    savefig('docs/images/noisy_sine_new_domain.png')


    figure(400)
    hold(True)
    plot(fit.x, fit.dy, '.-', color=belize_hole, linewidth=2)
    grid(True)
    xlabel('Time (seconds)')
    ylabel('Derivative of Amplitude (volts/sec)')
    ylim([-15, 15])
    savefig('docs/images/noisy_sine_fit_deriv.png')


