#!/usr/bin/env python
import numpy as np
import cmath as cm

class s_params:
    def __init__(self, s11, s12, s21, s22):
        self.s11 = s11
        self.s12 = s12
        self.s21 = s21
        self.s22 = s22

        self._gamma_out = None
        self._gamma_in  = None
        self._ucstable  = False

    def _det(self):
        return self.s11*self.s22 - self.s12*self.s21

    def rollet(self):
        det = self._det()
        num = 1 - np.abs(self.s11)**2 - np.abs(self.s22)**2 + np.abs(det)**2
        den = 2 * np.abs(self.s12*self.s21)
        k   = num/den
        if k > 1 and np.abs(det) < 1:
            self._ucstable = True
        return k, det

    def mu(self):
        '''Edwards-Sinsky single-parameter stability. Unconditionally stable if mu1 > 1 or mu2 > 1.'''
        det = self._det()
        mu1 = (1 - np.abs(self.s11)**2) / (np.abs(self.s22 - det*np.conj(self.s11)) + np.abs(self.s12*self.s21))
        mu2 = (1 - np.abs(self.s22)**2) / (np.abs(self.s11 - det*np.conj(self.s22)) + np.abs(self.s12*self.s21))
        return mu1, mu2

    def stability_circles(self):
        '''
        Returns output and input stability circle parameters.
          C_L, r_L : center and radius in Gamma_L (Smith chart) plane
          C_S, r_S : center and radius in Gamma_S plane
        Stable side: if |S11| < 1, the Smith chart origin is in the stable region for Gamma_L.
        '''
        det = self._det()

        C_L = np.conj(self.s22 - det*np.conj(self.s11)) / (np.abs(self.s22)**2 - np.abs(det)**2)
        r_L = np.abs(self.s12*self.s21) / np.abs(np.abs(self.s22)**2 - np.abs(det)**2)

        C_S = np.conj(self.s11 - det*np.conj(self.s22)) / (np.abs(self.s11)**2 - np.abs(det)**2)
        r_S = np.abs(self.s12*self.s21) / np.abs(np.abs(self.s11)**2 - np.abs(det)**2)

        return C_L, r_L, C_S, r_S

    def conjugate_match(self):
        '''
        Simultaneous conjugate match for maximum gain (valid only when unconditionally stable).
        Returns (Gamma_S_opt, Gamma_L_opt).
        '''
        det = self._det()

        B1 = 1 + np.abs(self.s11)**2 - np.abs(self.s22)**2 - np.abs(det)**2
        C1 = self.s11 - det*np.conj(self.s22)
        disc1 = np.sqrt(B1**2 - 4*np.abs(C1)**2 + 0j)
        gamma_s = (B1 - disc1) / (2*C1)
        if np.abs(gamma_s) > 1:
            gamma_s = (B1 + disc1) / (2*C1)

        B2 = 1 + np.abs(self.s22)**2 - np.abs(self.s11)**2 - np.abs(det)**2
        C2 = self.s22 - det*np.conj(self.s11)
        disc2 = np.sqrt(B2**2 - 4*np.abs(C2)**2 + 0j)
        gamma_l = (B2 - disc2) / (2*C2)
        if np.abs(gamma_l) > 1:
            gamma_l = (B2 + disc2) / (2*C2)

        return gamma_s, gamma_l

    def mag(self):
        '''
        Maximum available gain (MAG) if K >= 1, otherwise maximum stable gain (MSG).
        MAG = |S21/S12| * (K - sqrt(K^2 - 1))
        MSG = |S21/S12|
        '''
        k, _ = self.rollet()
        msg = np.abs(self.s21 / self.s12)
        if k >= 1:
            return msg * (k - np.sqrt(k**2 - 1))
        return msg

    def gmsg(self):
        return np.abs(self.s21/self.s12)

    def gamma_in(self, gamma_l):
        if self.s12 == 0 or self.s21 == 0:
            return self.s11
        num = self.s12*self.s21*gamma_l
        den = 1 - self.s22*gamma_l
        self._gamma_in = self.s11 + num/den
        return self._gamma_in

    def gamma_out(self, gamma_s):
        if self.s12 == 0 or self.s21 == 0:
            return self.s22
        num = self.s12*self.s21*gamma_s
        den = 1 - self.s11*gamma_s
        self._gamma_out = self.s22 + num/den
        return self._gamma_out

    def z_in(self, gamma_l, z0=50):
        gin = self.gamma_in(gamma_l)
        return z0 * (1 + gin) / (1 - gin)

    def z_out(self, gamma_s, z0=50):
        gout = self.gamma_out(gamma_s)
        return z0 * (1 + gout) / (1 - gout)

    def gain_t(self, gamma_s, gamma_l, gamma_out=None, gamma_in=None):
        gamma_out = gamma_out if gamma_out is not None else self._gamma_out
        gamma_in  = gamma_in  if gamma_in  is not None else self._gamma_in

        if self.s21 == 0:
            return 0
        if gamma_out is None and gamma_in is None:
            print("gain_t: call gamma_out(gamma_s) or gamma_in(gamma_l) first.")
            return -1
        if gamma_out is not None:
            lhs = (1 - np.abs(gamma_s)**2) / np.abs(1 - self.s11*gamma_s)**2
            rhs = (1 - np.abs(gamma_l)**2) / np.abs(1 - gamma_out*gamma_l)**2
            return lhs * np.abs(self.s21)**2 * rhs
        lhs = (1 - np.abs(gamma_s)**2) / np.abs(1 - gamma_in*gamma_s)**2
        rhs = (1 - np.abs(gamma_l)**2) / np.abs(1 - self.s22*gamma_l)**2
        return lhs * np.abs(self.s21)**2 * rhs

    def gain_p(self, gamma_l, gamma_in=None):
        gamma_in = gamma_in if gamma_in is not None else self._gamma_in

        if self.s21 == 0:
            return 0
        if gamma_in is None:
            print("gain_p: call gamma_in(gamma_l) first.")
            return -1
        lhs = 1 / (1 - np.abs(gamma_in)**2)
        rhs = (1 - np.abs(gamma_l)**2) / np.abs(1 - self.s22*gamma_l)**2
        return lhs * np.abs(self.s21)**2 * rhs

    def gain_a(self, gamma_s, gamma_out=None):
        gamma_out = gamma_out if gamma_out is not None else self._gamma_out

        if self.s21 == 0:
            return 0
        if gamma_out is None:
            print("gain_a: call gamma_out(gamma_s) first.")
            return -1
        lhs = (1 - np.abs(gamma_s)**2) / np.abs(1 - self.s11*gamma_s)**2
        rhs = 1 / (1 - np.abs(gamma_out)**2)
        return lhs * np.abs(self.s21)**2 * rhs
