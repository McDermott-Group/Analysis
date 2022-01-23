import numpy as np
import matplotlib.pyplot as plt

w = [1]
file = ["Q1(ri=90, ro=182).txt","Q2(ri=50, ro=52.9).txt","Q3_Blackboxcalc(with leads).txt","Q4(ri=70, ro=90.5).txt"]
q = 3

f = np.loadtxt("testpad_junc.txt", usecols=[0], skiprows=3)
f = f*1e9

Z_img = np.zeros([len(f),len(w)], dtype=float)
Z_real = np.zeros([len(f),len(w)], dtype=float)

for i in range(len(w)):
     Z_real[:,i] = np.loadtxt("testpad_junc.txt", usecols=[1], skiprows=3)
     Z_img[:,i] = np.loadtxt("testpad_junc.txt", usecols=[2], skiprows=3)


R = 5000
C = 3.5e-15
om = 2*np.pi*f
Z_JJ = np.zeros(len(f), dtype=complex)
Z_JJ = (1/R - 1j*om*C)/((1/R)**2 + (om*C)**2)
Z_rad = Z_real + 1j*Z_img
gamma = np.zeros([len(f),len(w)], dtype=complex)

for i in range(len(f)):
    gamma[i,:] = (Z_rad[i] - np.conj(Z_JJ[i]))/(Z_rad[i] + Z_JJ[i])

def complex_mod(z):
    a = z.real
    b = z.imag
    return np.sqrt(a**2 + b**2)

ec = 1 - (complex_mod(gamma))**2
ec_dB = 10*np.log10(ec)

plt.figure()
plt.plot(f, ec)
plt.title("Coupling effeciency_testpad (R={0} Ohms, C={1} F)".format(np.round(R),C))
plt.xlabel("frequency (GHz)")
plt.ylabel(r"coupling efficiency $e_c$ (dB)")
plt.show()



fa = np.loadtxt(file[q-1], usecols=[0], skiprows=3)
fa = fa*1e9

Z_imga = np.zeros([len(f),len(w)], dtype=float)
Z_reala = np.zeros([len(f),len(w)], dtype=float)

for i in range(len(w)):
     Z_reala[:,i] = np.loadtxt(file[q-1], usecols=[1], skiprows=3)
     Z_imga[:,i] = np.loadtxt(file[q-1], usecols=[2], skiprows=3)

if q != 3:
    Ra = 10000
    Ca = 3.5e-15
    oma = 2*np.pi*fa
    Z_JJa = np.zeros(len(f), dtype=complex)
    Z_JJa = (1/Ra - 1j*oma*Ca)/((1/Ra)**2 + (oma*Ca)**2)
    Z_rada = Z_reala + 1j*Z_imga
    gammaa = np.zeros([len(f),len(w)], dtype=complex)

    for i in range(len(f)):
        gammaa[i,:] = (Z_rada[i] - np.conj(Z_JJa[i]))/(Z_rada[i] + Z_JJa[i])

    def complex_mod(z):
        a = z.real
        b = z.imag
        return np.sqrt(a**2 + b**2)

    eca = 1 - (complex_mod(gammaa))**2
    ec_dBa = 10*np.log10(eca)

    plt.figure()
    plt.plot(fa, ec_dBa)
    plt.title("Coupling effeciency_Q{2} (Ra={0} Ohms, Ca={1} F)".format(np.round(Ra),Ca,q))
    plt.xlabel("frequency (GHz)")
    plt.ylabel(r"coupling efficiency $e_c$ (dB)")

else:
    Z_at = Z_reala + 1j*Z_imga

    Y_at = 1/Z_at


    Cj = 3.5e-15
    Lj = 21e-9
    Y_j = 1j*(2*np.pi*f)*Cj - 1j/(2*np.pi*f*Lj)

    Y_t = Y_at + Y_j

    Z_t = 1/Y_t

    plt.figure()
    plt.plot(f, Y_t.real)
    plt.plot(f, Y_t.imag, linestyle="dashed")
    plt.grid()


ec_t = ec*eca
ec_t_dB = 10*np.log10(ec_t)
plt.figure()
plt.plot(fa, ec_t_dB, label="ec_a * ec_r")
plt.title("Coupling effeciency_Q{0}".format(q))
plt.xlabel("frequency (GHz)")
plt.ylabel(r"coupling efficiency $e_c$ (dB)")
plt.legend()








































# C = [0.5e-15,2e-15,5e-15,10e-15]
# plt.figure()
# for j in range(len(C)):
#     R = 10000
#     C1 = C[j]
#     om = 2*np.pi*f
#     Z_JJ = np.zeros(len(f), dtype=complex)
#     Z_JJ = (1/R - 1j*om*C1)/((1/R)**2 + (om*C1)**2)
#     Z_rad = Z_real + 1j*Z_img
#     gamma = np.zeros([len(f),len(w)], dtype=complex)

#     for i in range(len(f)):
#         gamma[i,:] = (Z_rad[i] - np.conj(Z_JJ[i]))/(Z_rad[i] + Z_JJ[i])

#     def complex_mod(z):
#         a = z.real
#         b = z.imag
#         return np.sqrt(a**2 + b**2)

#     ec = 1 - (complex_mod(gamma))**2
#     ec_dB = 10*np.log10(ec)


#     plt.plot(f, ec_dB, label="C = {0}".format(C1))
#     plt.title("Coupling effeciency_Q4 (R={0} Ohms, C={1} F)".format(np.round(R),C))
#     plt.xlabel("frequency (GHz)")
#     plt.ylabel(r"coupling efficiency $e_c$ (dB)")
#     plt.legend()

















