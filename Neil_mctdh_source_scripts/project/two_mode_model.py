# system imports
# import os
# from os.path import join
import itertools as it
# import warnings

# third party imports
import scipy
from scipy.fftpack import fft
import numpy as np
# from numpy import newaxis as NEW
# import matplotlib as mpl
import matplotlib.pyplot as plt

# local imports
# from .vibronic import vIO, VMK


# factoring out the debug statements
print_2nd_quantized_hamiltonian = False
print_during_initialization = False
print_during_sos_routine = False


# just for testing / benchmarking
# one surface
class model_two_mode(object):
    "define a class that modeling "
    "dynamics of a harmonic potential"
    "with two vibrational modes"

    def __init__(self, E_0, L_Freq, omega_Ij, omega_ij, E_x=0., f=np.array([0., 0.])):
        "define a input model Hamiltonian with several hyperparameters"
        "set hbar=1, m=1"

        # set ground state energy
        self.f = f
        self.f_bar = np.ones(2) + f
        self.h_0 = E_0
        self.omega = omega_Ij.copy()
        self.omega_ij = omega_ij.copy()
        self.omega_IJ = omega_ij.copy()
        self.omega_i = L_Freq.copy()
        self.omega_I = L_Freq.copy()

        # diagnostic printing?
        if print_2nd_quantized_hamiltonian:
            print(
                '***2nd quantitized model Hamiltonian***',
                '***h_0', self.h_0, '\n',
                '***omega', self.omega, '\n',
                '***omega_i', self.omega_i, '\n',
                '***omega_I', self.omega_I, '\n',
                '***omega_ij', self.omega_ij, '\n',
                '***omega_IJ', self.omega_IJ, '\n',
                '*** End of Initialization *** ',
                sep="\n"
            )

        # define residue equation for I
        def f_0(t, t_I, t_i, t_IJ, t_ij):
            R = self.h_0 + 0.0*1j  # why?
            for i, j in it.product(range(2), range(2)):
                R += self.omega[i, j] * t[j, i]
                R += self.omega[i, j] * t_I[i] * t_i[j]
                R += 0.5 * self.omega_ij[i, j] * t_IJ[i, j]
                R += 0.5 * self.omega_ij[i, j] * t_I[i] * t_I[j]
                R += 0.5 * self.omega_IJ[i, j] * t_ij[i, j]
                R += 0.5 * self.omega_IJ[i, j] * t_i[i] * t_i[j]
            for i in range(2):
                R += self.omega_i[i] * t_I[i]
                R += self.omega_I[i] * t_i[i]
            return R

        self.f_0 = f_0

        # define residue equation for i^dagger
        def f_1(t, t_I, t_i, t_IJ, t_ij):
            R = np.zeros(2, dtype=complex)
            for e in range(2):
                for i, j in it.product(range(2), range(2)):
                    R[e] += self.omega[i, j] * t_I[i] * t[j, e]
                    R[e] += self.omega[i, j] * t_i[j] * t_IJ[i, e]
                    R[e] += self.omega_ij[i, j] * t_I[i] * t_IJ[j, e]
                    R[e] += self.omega_IJ[i, j] * t_i[j] * t[i, e]
                for i in range(2):
                    R[e] += self.omega_i[i] * t_IJ[i, e]
                    R[e] += self.omega_I[i] * t[i, e]
                    R[e] += self.f[e] * self.omega[i, e] * t_I[i]
                    R[e] += self.f[e] * self.omega_IJ[e, i] * t_i[i]
                R[e] += self.f[e] * self.omega_I[e]
            return R
        self.f_1 = f_1

        # define residue equation for i
        def f_2(t, t_I, t_i, t_IJ, t_ij):
            R = np.zeros(2, dtype=complex)
            for e in range(2):
                for i, j in it.product(range(2), range(2)):
                    R[e] += self.omega[i, j] * t_i[j] * t[e, i]
                    R[e] += self.omega[i, j] * t_I[i] * t_ij[j, e]
                    R[e] += self.omega_ij[i, j] * t_I[i] * t[e, j]
                    R[e] += self.omega_IJ[i, j] * t_i[i] * t_ij[j, e]
                for i in range(2):
                    R[e] += self.f_bar[e] * self.omega[e, i] * t_i[i]
                    R[e] += self.omega_i[i] * t[e, i]
                    R[e] += self.omega_I[i] * t_ij[i, e]
                    R[e] += self.f_bar[e] * self.omega_ij[i, e] * t_I[i]
                R[e] += self.f_bar[e] * self.omega_i[e]
            return R
        self.f_2 = f_2

        # define residue equation for i^daggerj
        def f_3(t, t_I, t_i, t_IJ, t_ij):
            R = np.zeros((2, 2), dtype=complex)
            for e, f in it.product(range(2), range(2)):
                for i, j in it.product(range(2), range(2)):
                    R[e, f] += self.omega[i, j] * t[f, i] * t[j, e]
                    R[e, f] += self.omega_IJ[i, j] * t_ij[i, f] * t[j, e]
                    R[e, f] += self.omega_ij[i, j] * t_IJ[i, e] * t[f, j]
                for i in range(2):
                    R[e, f] += self.f_bar[f] * self.omega_ij[i, f] * t_IJ[i, e]
                    R[e, f] += self.f_bar[f] * self.omega[f, i] * t[i, e]
                    R[e, f] += self.f[e] * self.omega_IJ[e, i] * t_ij[i, f]
                R[e, f] += self.f[e] * self.f_bar[f] * self.omega[f, e]
            return R
        self.f_3 = f_3

        # define residue equation for i^dagger j^dagger
        def f_4(t, t_I, t_i, t_IJ, t_ij):
            R = np.zeros([2, 2], dtype=complex)
            for e, f in it.product(range(2), range(2)):
                for i, j in it.product(range(2), range(2)):
                    R[e, f] += self.omega[i, j] * t[j, f] * t_IJ[i, e]
                    R[e, f] += self.omega[i, j] * t[j, e] * t_IJ[i, f]
                    R[e, f] += self.omega_ij[i, j] * t_IJ[j, f] * t_IJ[i, e]
                    R[e, f] += self.omega_IJ[i, j] * t[i, e] * t[j, f]
                for i in range(2):
                    R[e, f] += self.f[f] * self.omega[i, f] * t_IJ[i, e]
                    R[e, f] += self.f[f] * self.omega_IJ[f, i] * t[i, e]
                    R[e, f] += self.f[e] * self.omega_IJ[e, i] * t[i, f]
                    R[e, f] += self.f[e] * self.omega[e, i] * t[i, f]
                R[e, f] += self.f[e] * self.f[f] * self.omega_IJ[e, f]
            return R
        self.f_4 = f_4

        # define residue equation for ij
        def f_5(t, t_I, t_i, t_IJ, t_ij):
            R = np.zeros([2, 2], dtype=complex)
            for e, f in it.product(range(2), range(2)):
                for i, j in it.product(range(2), range(2)):
                    R[e, f] += self.omega_ij[i, j] * t[f, i] * t[e, j]
                    R[e, f] += self.omega_IJ[i, j] * t_ij[j, e] * t_ij[i, f]
                    R[e, f] += self.omega[i, j] * t_ij[j, e] * t[f, i]
                    R[e, f] += self.omega[i, j] * t_ij[j, f] * t[e, i]
                for i in range(2):
                    R[e, f] += self.f_bar[f] * self.omega[f, i] * t_ij[i, e]
                    R[e, f] += self.f_bar[e] * self.omega[e, i] * t_ij[i, f]
                    R[e, f] += self.f_bar[f] * self.omega_ij[f, i] * t[e, i]
                    R[e, f] += self.f_bar[e] * self.omega_ij[e, i] * t[f, i]
                R[e, f] += self.f_bar[e] * self.f_bar[f] * self.omega_ij[e, f]
            return R
        self.f_5 = f_5

    def solve(self, t_init=0., t_term=10., N=300000, f_=np.zeros(2), E_g=0., BF=1.):
        "solve solve numerical differential of the two mode systems"
        self.f = f_
        self.f_bar = np.ones(2) + f_
        # self.h_0-=self.f[0]*self.omega[0,0]+self.f[1]*self.omega[1,1]
        print("***f= " + str(self.f))
        print("***f_bar= " + str(self.f_bar))
        # if self.f[0]==0 and self.f[1]==0:
        #     self.omega[0,0]-=self.omega_g1
        #     self.omega[1,1]-=self.omega_g2
        # self.h_0+=0.5*(self.omega_g1+self.omega_g2)
        # self.omega_IJ=2*self.omega_IJ
        # self.omega_ij=2*self.omega_ij
        step = (t_term - t_init) / N
        self.T_step = np.linspace(t_init, t_term, N)
        self.T_step = self.T_step
        C = 0.0 + 0.0*1j
        t_I = np.zeros(2, dtype=complex)
        t_i = np.zeros(2, dtype=complex)
        t = np.zeros([2, 2], dtype=complex)
        t_IJ = np.zeros([2, 2], dtype=complex)
        t_ij = np.zeros([2, 2], dtype=complex)
        t_I_step = []
        t_i_step = []
        t_step = []
        t_IJ_step = []
        t_ij_step = []
        C_step = []
        # print(t_I-self.f_1(t,t_I,t_i,t_IJ,t_ij)*step*1j)
        for _ in range(N):
            # update auto-correlation function
            # correct=0
            # for i in range(2):
            #     correct+=self.f[i]*t[i,i]
            C -= (self.f_0(t, t_I, t_i, t_IJ, t_ij) + self.E_x) * step * 1j  # -correct
            # update amplitude
            t_I_old = t_I.copy()
            t_i_old = t_i.copy()
            t_old = t.copy()
            t_ij_old = t_ij.copy()
            t_IJ_old = t_IJ.copy()
            t_I -= self.f_1(t_old, t_I_old, t_i_old, t_IJ_old, t_ij_old) * step * 1j
            t_i -= self.f_2(t_old, t_I_old, t_i_old, t_IJ_old, t_ij_old) * step * 1j
            t -= self.f_3(t_old, t_I_old, t_i_old, t_IJ_old, t_ij_old) * step * 1j
            t_IJ -= self.f_4(t_old, t_I_old, t_i_old, t_IJ_old, t_ij_old) * step * 1j
            t_ij -= self.f_5(t_old, t_I_old, t_i_old, t_IJ_old, t_ij_old) * step * 1j

            # update auto-correlation function
            # C-=self.f_0(t,t_I,t_i,t_IJ,t_ij)*step*1j

            # store amplitude
            t_I_step.append(t_I)
            t_i_step.append(t_i)
            t_step.append(t)
            t_IJ_step.append(t_IJ)
            t_ij_step.append(t_ij)

            # store energy
            C_step.append(C)
            if _ % 10000 == 0:
                print(C_step[_])
        self.t_I = t_I_step
        self.t_i = t_i_step
        self.t = t_step
        self.t_IJ = t_IJ_step
        self.t_ij = t_ij_step
        self.C_step = np.exp(np.array(C_step))

        # Fourier transform the time correlation function
        I_omega = np.fft.fft(self.C_step) / self.C_step.size
        self.I_omega = I_omega * BF
        self.freq = 2. * np.pi * np.fft.fftfreq(self.C_step.size, step)
        # print(T)
        # self.F_range = np.linspace(0.,1./(T),100000)
        return self.T_step, self.C_step

    def sos_solution(self, HO_size=40, flag=[0, 0], E_g=0, BF=1):
        "Solve the two mode H.O. problem by Sum Over States (SOS) scheme"
        # define the matrix element of Hamiltonian in H.O. basis
        # E_g=self.E_g[0]
        # flag=[0,0]

        if flag[0] == 0 and flag[1] == 0:
            self.omega_ij = 0.5 * self.omega_ij
            self.omega_IJ = 0.5 * self.omega_IJ
            # self.omega[0,0] += self.omega_g1
            # self.omega[1,1] += self.omega_g2
            # self.h_0 += 0.5*(self.omega_g1 + self.omega_g2)
            # self.h_0 -= (self.f[0]*self.omega[0,0] + self.f[1]*self.omega[1,1])

        H = np.zeros([HO_size, HO_size, HO_size, HO_size], dtype=complex)

        for a_1, a_2, b_1, b_2 in it.product(
            range(HO_size), range(HO_size), range(HO_size), range(HO_size)
        ):

            index = (a_1, a_2, b_1, b_2)

            if a_1 == b_1 and a_2 == b_2:
                H[index] = self.h_0
                H[index] += self.omega[0, 0] * (b_1)
                H[index] += self.omega[1, 1] * (b_2)

            if a_1 == b_1+1 and a_2 == b_2-1:
                H[index] = self.omega[0, 1] * np.sqrt(b_1+1) * np.sqrt(b_2)

            if a_1 == b_1-1 and a_2 == b_2+1:
                H[index] = self.omega[1, 0] * np.sqrt(b_1) * np.sqrt(b_2+1)

            if a_1 == b_1+1 and a_2 == b_2:
                H[index] = self.omega_I[0] * np.sqrt(b_1+1)

            if a_1 == b_1 and a_2 == b_2+1:
                H[index] = self.omega_I[1] * np.sqrt(b_2+1)

            if a_1 == b_1-1 and a_2 == b_2:
                H[index] = self.omega_i[0] * np.sqrt(b_1)

            if a_1 == b_1 and a_2 == b_2-1:
                H[index] = self.omega_i[1] * np.sqrt(b_2)

            if a_1 == b_1+2 and a_2 == b_2:
                H[index] = self.omega_IJ[0, 0] * np.sqrt(b_1+1) * np.sqrt(b_1+2)

            if a_1 == b_1+1 and a_2 == b_2+1:
                H[index] = (self.omega_IJ[0, 1] + self.omega_IJ[1, 0]) * \
                    np.sqrt(b_1+1) * np.sqrt(b_2+1)

            if a_1 == b_1 and a_2 == b_2+2:
                H[index] = self.omega_IJ[1, 1] * np.sqrt(b_2+1) * np.sqrt(b_2+2)

            if a_1 == b_1-2 and a_2 == b_2:
                H[index] = self.omega_ij[0, 0] * np.sqrt(b_1) * np.sqrt(b_1-1)

            if a_1 == b_1-1 and a_2 == b_2-1:
                H[index] = (self.omega_ij[0, 1] + self.omega_ij[1, 0]) * \
                    np.sqrt(b_1) * np.sqrt(b_2)

            if a_1 == b_1 and a_2 == b_2-2:
                H[index] = self.omega_ij[1, 1] * np.sqrt(b_2) * np.sqrt(b_2-1)

        # reshape the tensor to a square matrix
        H = H.reshape(HO_size**2, HO_size**2)

        if print_during_sos_routine:
            print(
                "check Hermicity:",
                abs(H-H.transpose()).max(),
                "matrix element of Hamiltonian",
                H,
                sep="\n",
            )

        # always make sure each block is Hermitian
        assert np.allclose(abs(H - H.transpose()), 0.0), "This block isn't Hermitian"

        # diagonalize the Hamiltonian
        E, V = np.linalg.eigh(H)

        if print_during_sos_routine:
            print("Eigen Energies")
            for i in range(10):
                print(E[i])

        return H

        # test eigen energies:
        # E_test=[]
        # l,C=np.linalg.eigh(self.W)  # construct vibrational mode in orthornormal basis
        # omega=np.sqrt(l)
        # print("orthornormal vibrational modes:")
        # print(omega)
        # for n_1 in range(0,40):
        #     for n_2 in range(0,40):
        #         E_test.append(omega[0]*(n_1+0.5)+omega[1]*(n_2+0.5)-sum(self.g**2/(2*l)))
        # E_test=np.sort(np.array(E_test))
        # print("expected eigen energies")
        # for i in range(10):
        #     print(E_test[i]-E_g)

        # construct time correlation function

        # t=np.linspace(0.,10.,100000)
        # C_tau=np.zeros_like(t,dtype=complex)
        # norm=0
        # band=40*flag[0]+flag[1]
        # for n in range(10):
        #     norm+=V[band,n]**2
        #     C_tau+=np.exp(-1j*(E[n]+self.E_x)*t)*V[band,n]**2
        # self.E_zero=E.min() #get zero point energy
        # C_tau=C_tau/norm
        # self.C_tau_sos=C_tau
        # self.t_sos=t
        # plt.plot(t,C_tau.real,label='real')
        # plt.plot(t,C_tau.imag,label='imag')
        # plt.title("Plot of time correlation function C(t) vs t")
        # plt.xlabel("t")
        # plt.ylabel("C(t)")
        # plt.ylim(-100,100)
        # plt.legend()
        # plt.show()
        # yf = fft(C_tau)/C_tau.size
        # freq = 2*np.pi*np.fft.fftfreq(C_tau.size,d=0.0001)
        # print(freq)
        # self.freq_sos=freq
        # self.yf=np.abs(yf)*BF
        # return t,freq,self.C_tau_sos,self.yf
        # return t_tau,freq,C_tau,yf
        # plt.plot(freq, np.abs(yf)/C_tau.size)
        # plt.ylim(0,5e-2)
        # plt.xlim(-100,100)
        # plt.title("Power Spectrum of Frank Cordon Model Hamiltonian",fontsize=20)
        # plt.xlabel(r'$\omega$',fontsize=18)
        # plt.ylabel(r'$I(\omega)$',fontsize=18)
        # plt.show()