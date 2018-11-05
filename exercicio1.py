# 1) Implementar núcleo de da transformada de Fourier (2.5 Pontos)
# a) Versão DFT
# b) Versão FFT
# c) Executar a transformada para um sinal simulado de senoidais. Em 30Hz,50Hz e 60Hz

import sys
import numpy as np

def main(option):
    if option == 'DFT':
        DFT()
    elif option == 'FFT':
        FFT() 
    elif option == 'TRANSFORMADA':
        pass

if __name__ == '__main__':

    main(sys.argv[1])

def DFT(fnlist):
    N = len(fnList)
    FmList = []
    for m in range(N):
        Fm = 0.0
        for n in range(N):
            Fm += fnList[n] * cmath.exp(- 1j * pi2 * m * n / N)
        FmList.append(Fm / N)
    return FmList

def FFT(x):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    
    if N % 2 > 0:
        raise ValueError("size of x must be a power of 2")
    elif N <= 32:  # this cutoff should be optimized
        return DFT_slow(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:N / 2] * X_odd,
                               X_even + factor[N / 2:] * X_odd])