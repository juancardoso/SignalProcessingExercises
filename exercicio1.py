# 1) Implementar núcleo de da transformada de Fourier (2.5 Pontos)
# a) Versão DFT
# b) Versão FFT
# c) Executar a transformada para um sinal simulado de senoidais. Em 30Hz,50Hz e 60Hz
# 2) Plote o gráfico do espectro dos dados colhidos (2.5 Pontos)
# 3) Usando o microfone faça um processo de aquisição do som em “tempo real” e plote o espectro a
# cada 200ms. (2.5 Pontos). Dica utilize slimDX
# 4) Com os dados do item 3, faça o desenho os espectrograma no tempo. Também plotando um conjunto
# a cada 200ms. (2.5 Pontos). Dica: Utilize o slimDX.


import sys
import numpy as np
import math
import random
import cmath
import matplotlib.pyplot as plt


Fs = 128.0
Ts = 1.0/Fs
t = np.arange(0,1,Ts)
pi2 = cmath.pi * 2.0

def loadList(N,f):
    a = float(random.randint(1, 100))
    fnList = []
    for n in range(N):
        t = float(n) / N * pi2
        fn = a * math.sin(f * t + float(random.randint(0, 360) / 360 * pi2))
        fnList.append(fn)
    return fnList

def DFT(fnList):
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
        return DFT(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:int(N / 2)] * X_odd,
                               X_even + factor[int(N / 2):] * X_odd])

def plot(x, fnList):
    fn = lambda x : math.sqrt(math.pow(x[0], 2) + math.pow(x[1], 2))
    X = fn(fnList)
    Y = fnList[range(int(len(fnList)/2))]

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(t,x)
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('Amplitude')
    ax[1].plot(frq,abs(Y),'r') # plotting the spectrum
    ax[1].set_xlabel('Freq (Hz)')
    ax[1].set_ylabel('|Y(freq)|')
    plt.show()

def main(option):
    if option == 'DFT':
        DFT()
    elif option == 'FFT':
        FFT() 
    elif option == 'TRANSFORMADA':
        n = int(sys.argv[2])
        f = int(sys.argv[3])
        lst = loadList(n, Fs)
        fft = FFT(lst)
        print(fft)
        plot(lst, fft)

if __name__ == '__main__':

    main(sys.argv[1])