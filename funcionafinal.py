import numpy as np
from numpy.fft import fft
import scipy.signal as signal
import scipy.linalg as linalg
import scipy.io.wavfile as wavfile
from sympy import *
from IPython.display import display
import scipy.signal as signal
import matplotlib.pyplot as plt

# escritura del vector
fs, x = wavfile.read('C:\\Users\\santo\\Desktop\\TAREAFINAL COM4\\x3_U017.wav')
lx = len(x)
fs, y = wavfile.read('C:\\Users\\santo\\Desktop\\TAREAFINAL COM4\\y3_U017.wav')
ly = len(y)

# creacion del vector para graficar señal de audio
n = np.arange(lx)
import matplotlib.pyplot as plt
plt.subplot(2,1,1)
plt.plot(n,x)
plt.title('X[n]')
plt.subplot(2,1,2)
plt.stem(n,y)
plt.title('Y[n]')
plt.show()

# para solucion de problema
m = int(input('Ingrese el valor limite para los coeficientes b (M): '))
k = int(input('Ingrese el valor limite para los coeficientes a (N): '))

# Determinación de numero de incognitas
M = m+2
m = m+2
k = m+k-1
q = np.arange(3000) # limite para que lo solucione más rápido
if k > q.size:
    print('ERROR el valor de N y M es mayor que la longitud del vector')
else:
    # creacion de coeficientes del sistema de ecuaciones
    C = np.zeros((q.size, k))
    bl = np.zeros(q.size)

    for j in range(q.size):
        l = 1
        for i in range(k):
            if i < M:
                C[j, i] = x[m-i]
            if i > M-1:
                C[j, i] = -y[m-l-1]
                l = l + 1
        bl[j] = y[m-1]
        
        m = m + 1

    B = np.transpose(bl)
    R = linalg.lstsq(C, B)[0] # Solucion a traves de minimos cuadrados

    #R = linalg.lstsq(C, B, rcond=None)[0] # Solucion a traves de minimos cuadrados
    print('EL valor para cada uno de los coeficientes es: [b_0; b_1 .......;b_M; a_o;..... A_N]')
    print(R)

R = np.transpose(R)

# Arreglo de vectores para graficos
M = int(M)
b = R[0:M-1]
a = np.zeros(k - M + 2)
for j in range(M, k+1):
    a[j-M+1] = R[j-M+2]


# Encontrando el valor de H
z = symbols('z', real=True)

[r, p, k] = signal.residue(b, a)
HH = 0
for i in range(q.size):
    freq = q[i]/len(q) # normalización de la frecuencia
    #HH += sum(r/(1-p*freq)) + sum(k*np.array([1, 1/freq]))
    HH += sum(r/(1-p*freq)) + sum(k*np.array([1, 1/freq])) if freq > 1e-12 else 0


print('La función H (con la fracción parciales) para el sistema es:')
display(HH)

# ANALISIS ESPECTRAL
f = 0.1 # frecuencia de 0.1 Hz
N1 = 30
N2 = 120

# Analisis Espectral para x
# transformadas:
X1 = abs(np.fft.fft(x, N1))
X2 = abs(np.fft.fft(x, N2))

# Rango normalizado para transformadas:
F1x = np.arange(N1)/N1
F2x = np.arange(N2)/N2

# Analisis Espectral para Y
Y1 = abs(np.fft.fft(y, N1))
Y2 = abs(np.fft.fft(y, N2))

# Rango normalizado para transformadas:
F1y = np.arange(N1)/N1
F2y = np.arange(N2)/N2

# Grafica de funciones
# Funciones X y Y
plt.subplot(4, 2, 1)
plt.stem(n, x)
plt.title('X[n]')

plt.subplot(4, 2, 2)
plt.stem(n, y)
plt.title('Y[n]')

# Graficas del espectro:
plt.subplot(4, 2, 3)
plt.stem(F1x, X1, '.')
plt.title('Espectro par X[n] con N=30')

plt.subplot(4, 2, 5)
plt.stem(F2x, X2, '.')
plt.title('Espectro par X[n] con N=120')

plt.subplot(4, 2, 4)
plt.stem(F1y, Y1, '.')
plt.title('Espectro para Y[n] con N=30')

plt.subplot(4, 2, 6)
plt.stem(F2y, Y2, '.')
plt.title('Espectro para Y[n] con N=120')

#plt.show()
# Evaluar la función HH en un conjunto de puntos
z = symbols('z', real=True)
HH_eval = lambdify(z, HH, 'numpy')

HH_real = np.zeros(q.size)
HH_imag = np.zeros(q.size)
for i in range(q.size):
    HH_val = HH_eval(q[i])
    HH_real[i] = HH_val.real
    HH_imag[i] = HH_val.imag

plt.stem(q, HH_real)
plt.title('Parte real H[Z]')
plt.show()

plt.stem(q, HH_imag)
plt.title('Parte Imaginaria H[Z]')
plt.show()

# respuesta en amplitud y frecuencia
w, h = signal.freqz(b, a)
fig, ax1 = plt.subplots()
ax1.set_title('Respuesta en amplitud y fase de H[z]')
ax1.plot(w, 20*np.log10(abs(h)), 'b')
ax1.set_ylabel('Amplitud [dB]', color='b')
ax1.set_xlabel('Frecuencia [rad/muestra]')
ax1.grid()

ax2 = ax1.twinx()
angles = np.unwrap(np.angle(h))
ax2.plot(w, angles, 'g')
ax2.set_ylabel('Fase (radianes)', color='g')
ax2.grid()
plt.show()

# Polos y ceros
zeros = np.roots(b)
poles = np.roots(a)
fig, ax = plt.subplots()
ax.scatter(np.real(zeros), np.imag(zeros), marker='o', color='b')
ax.scatter(np.real(poles), np.imag(poles), marker='x', color='r')
plt.title('Diagrama de polos y cero')
plt.xlabel('Parte Real')
plt.ylabel('Parte Imaginaria')
plt.show()

import scipy.signal as signal
import matplotlib.pyplot as plt

# Obtener la respuesta al impulso
b, a = signal.butter(4, 0.2, 'lowpass')
h, t = signal.impulse((b, a))

# Graficar la respuesta al impulso
plt.stem(t, h)
plt.title('Respuesta al impulso')
plt.xlabel('n')
plt.ylabel('h[n]')
plt.show()

import scipy.signal as signal
import matplotlib.pyplot as plt

# Obtener la respuesta al escalón unitario
b, a = signal.butter(4, 0.2, 'lowpass')
r, t = signal.step((b, a))

# Graficar la respuesta al escalón unitario
plt.stem(t, r)
plt.title('Respuesta al escalón unitario')
plt.xlabel('n')
plt.ylabel('r[n]')
plt.show()

""" Este código en Python carga dos archivos de audio y los grafica en dos subplots diferentes. 
Luego, se le pide al usuario que ingrese dos valores para los límites de los coeficientes 
b y a, y se crea un sistema de ecuaciones que utiliza estos valores para encontrar la solución 
a través de mínimos cuadrados. Los coeficientes encontrados se utilizan para encontrar la función 
H del sistema y se grafica su parte real e imaginaria. También se grafica la respuesta en amplitud 
y fase de la función H, así como el diagrama de polos y ceros. Finalmente, se obtiene y grafica la 
respuesta al impulso y la respuesta al escalón unitario de un filtro de Butterworth de orden 4 y
frecuencia de corte de 0.2. """