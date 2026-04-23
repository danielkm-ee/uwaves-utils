#!/usr/bin/env python
# helper functions for conversions, etc
import numpy as np


def dbv(mag):
    return 20*np.log10(np.abs(mag))

def db(mag):
    if mag == 0:
        return -np.inf
    else:
        return 10*np.log10(np.abs(mag))

def polar(mag, deg):
    return mag * np.exp(1j * np.radians(deg))

def gamma2z(gamma, z0=50):
    return z0 * (1 + gamma) / (1 - gamma)

def z2gamma(z, z0=50):
    return (z - z0) / (z + z0)

def denormz(z, z0=50):
    return z*z0

def denormy(y, z0=50):
    return y/z0
