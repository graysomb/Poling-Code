import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def duty(x, Lcryst, cheiff, absDmax):
    # This function calculates the duty cycle to achieve a certain effective
    # nonlinearity. The duty cycle is limited to some absolute duty cycle
    # (-1,0,1) -> (0,100%)
    cheiffx = cheiff(x)
    d1 = -np.sqrt(1 - cheiffx)
    d1[d1 < -absDmax] = -absDmax
    d2 = np.sqrt(1 - cheiffx)
    d2[d2 > absDmax] = absDmax
    d = (x < 0) * d1 + (x > 0) * d2
    return d


if __name__ == "__main__":
    # This code calculates the phase matching function for a periodicaly poled KTP and a gaussian pump
    #currently it is for a .6 nm bandwidth gaussian pump at 775 nm undergoing degenerate SPDC.
    # The code generates a crystal using an analytic function, then calculates the widths for each
    #poling region. Finally the FFT is taken to determine the phase matching function.

    # Define units
    mm = 1e-3
    um = 1e-6
    nm = 1e-9
    c = 3e8

    # KTP refractive index from SPDcalc
    ktpx = lambda x: np.sqrt(2.10468 + 0.89342 * x**2 / (x**2 - 0.04438) - 0.01036 * x**2)
    ktpy = lambda x: np.sqrt(2.0993 + 0.922683 * x**2 / (x**2 - 0.0467695) - 0.0138408 * x**2)
    ktpz = lambda x: np.sqrt(1.9446 + 1.3617 * x**2 / (x**2 - 0.047) - 0.01491 * x**2)
    # Initialize pump parameters and functions
    lp = 0.775e-6  # pump wavelength
    wp = 2 * np.pi * c / lp
    kp = 2 * np.pi * ktpy(lp / um) / lp
    lam = lambda w: 2 * np.pi * c / w
    # convert to FWHM
    pbi = (2 / (2 * np.sqrt(2 * np.log(2)))) * 0.6e-9  # pump bandwidth
    pb = (2 * np.pi * c * pbi) / lp**2  # convert to freq
    aw = lambda ws, wi: np.exp(-((ws + wi - wp) / pb)**2)  # pump envelope function

    # Initialize crystal parameters and function
    k = lambda w, n: 2 * np.pi * n(lam(w) / um) / lam(w)
    dk = lambda ws, wi: k(ws + wi, ktpy) - k(ws, ktpy) - k(wi, ktpz)  # phase matching
    Lcryst = 27.38 * mm  # crystal length
    LamC = np.abs(2 * np.pi / dk(wp / 2, wp / 2))  # poling period

    # Create wavelength grid around degenerate 1550nm
    Nlam = 100  # wavelength points
    dLam = 10  # wavelength extent
    LamAll = np.linspace(1550 - dLam / 2, 1550 + dLam / 2, Nlam) * nm
    wall = 2 * np.pi * c / LamAll
    Ws, Wi = np.meshgrid(wall, np.flip(wall))

    # Create crystal in spatial domain and take FFT
    # FFT and sampling rate have to be multiples of each other
    N = 50  # samples per poling period
    LN = 20  # crystal lengths in entire spatial window (pads with zeros)
    Lall = Lcryst * LN
    Fs = 1 / (LamC / N)  # Sampling frequency
    T = 1 / Fs  # Sampling period
    L = int(Lall * Fs)  # Length of signal
    t = np.arange(0, L) * T  # Time vector
    polL = LamC  # poling period

    t = t - Lall / 2  # shift crystal position to zero

    # Calculate maximum duty cycle based on minimum poling period
    # mind = 0
    mind = 3.425e-6
    mind = mind / 4  # not sure why this is required
    absDmax = 1 - mind / (LamC)

    # Create apodization function
    FWHM = 14.60 * mm  # apodization width
    apd = FWHM / (2 * np.sqrt(2 * np.log(2)))  # convert to 1/e
    cheiff = lambda x: np.exp(-x**2 / apd**2)  # effective nonlinearity

    # Construct poling amplitude in spatial domain
    tN = lambda x: (x * Fs + Lcryst * Fs / 2)  # points in space
    Nw = int(Lcryst * Fs)  # number of points
    # function which outputs a square wave by taking the sign of a sin wave
    # with duty cycle determine by the effective nonlinearity
    sW = lambda x, f: np.sign(np.sin(x * (2 * np.pi * f)) + duty(x, Lcryst, cheiff, absDmax))

    # Calculate function over spatial domain
    S = sW(t, (1 / polL))
    # Cut out function over crystal region
    winC = np.logical_and(t > -Lcryst / 2, t < Lcryst / 2)

    # Convert crystal to a series of domain widths
    ti = np.where(np.diff(S * winC) != 0)[0] + 1
    widths = np.diff(ti) * T
    # Round the widths to .025um (given by raicol)
    widths = widths / (.025 * um)
    widths = np.round(widths) * .025 * um
    crys = np.array([])

    # Construct list of widths and convert back to the spatial domain
    for i in range(widths.size):
        sign = 1 - 2 * (i % 2)
        crys = np.hstack((crys, sign * np.ones(int(round(widths[i] / T)))))
        crys = np.hstack((crys, np.zeros(int(round(.025 * um / T)))))


    # Pad the crystal
    PadN = t.size - crys.size
    cryst = np.hstack((np.zeros(int(PadN / 2)), crys, np.zeros(int(PadN / 2))))

    plt.figure(2)
    plt.subplot(1, 2, 1)
    plt.plot(t, S * winC, t, cryst)  # plot both crystals, cheiff, and duty cycle

    plt.plot(t, duty(t, Lcryst, cheiff, absDmax))
    plt.plot(t, cheiff(t))
    plt.legend(['polling', 'duty', 'cheiff'])
    plt.xlim([-Lcryst, Lcryst])
    plt.title('Crystal Poling')
    plt.xlabel('z (m)')
    plt.ylabel('chi(z)')

    # Take FFT of crystal
    # Y = np.fft.fft(S * winC)
    Y = np.fft.fft(cryst)
    P2 = np.abs(Y / L)
    # Conversion due to DFT
    P1 = P2[:L // 2 + 1]
    P1[1:-1] = 2 * P1[1:-1]
    # Define spatial frequencies
    f = Fs * np.arange(0, L // 2 + 1) / L

    # Plot phase matching function
    plt.figure(2)
    plt.subplot(1, 2, 2)
    plt.plot(2 * np.pi * f, P1)

    plt.plot(2 * np.pi * f, np.max(P1) * np.abs(np.sinc(Lcryst * (2 * np.pi * f - 2 * np.pi / LamC) / (2 * np.pi))))
    plt.title('Phase Matching Function')
    plt.xlabel('dk (1/m)')
    plt.ylabel('h(dk)')
    plt.xlim([2 * np.pi / LamC * (1 - .01), 2 * np.pi / LamC * (1 + .01)])

    # Define an interpolation function for the phase matching curve
    phiwF = lambda ws, wi: interp1d(2 * np.pi * f, P1, fill_value='extrapolate')(-dk(ws, wi))

    # Calculate JSI
    JSA2 = phiwF(Ws, Wi) * aw(Ws, Wi)
    s2 = np.linalg.svd(JSA2, compute_uv=False)
    lambda2 = s2 ** 2
    lambda2 = lambda2 / np.sum(lambda2)
    K2 = 1 / np.sum(lambda2 ** 2)
    Leff = np.trapz(cheiff(t) * winC, t) / Lcryst

    # Normalize JSA and plot
    lambGs = 2 * np.pi * c / Ws - 1.55 * um
    lambGi = 2 * np.pi * c / Wi - 1.55 * um
    JSI = JSA2 ** 2 / np.max(np.max(JSA2 ** 2))
    JSI = JSA2 / np.max(np.max(JSA2))

    JSI = JSA2 / np.max(np.max(JSA2))
    f = plt.figure()
    plt.contourf(lambGs,lambGi,1 - JSI, origin='lower', cmap='bone')
    plt.xlabel('Signal (µm)')
    plt.ylabel('Idler (µm)')
    plt.title('JSA (Apod)')
    np.savetxt('poling_domain_widths.txt', widths/um)
    print("Schimdt Number: ", K2)
    plt.show()
