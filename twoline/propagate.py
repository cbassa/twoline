#!/usr/bin/env python3
import numpy as np
from twoline.twoline import TwoLineElement, set_checksum
from twoline.twoline import datetime_to_epoch, format_tle
from sgp4.api import Satrec, jday

def propagate(tle, t):
    # New epoch
    epochyr, epochdoy = datetime_to_epoch(t)
    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond / 1000000)
    
    # Compute reference state vector
    sat = Satrec.twoline2rv(tle.line1, tle.line2)
    e, r, v = sat.sgp4(jd, fr)
    r0, v0 = np.asarray(r), np.asarray(v)

    # Start loop
    rnew, vnew = r0, v0
    converged = False
    for i in range(100):
        # Convert state vector into new TLE
        newtle = classical_elements(rnew, vnew, epochyr, epochdoy, tle)
        
        # Propagate with SGP4
        sat = Satrec.twoline2rv(newtle.line1, newtle.line2)
        e, rtmp, vtmp = sat.sgp4(sat.jdsatepoch, sat.jdsatepochF)
        rsgp4, vsgp4 = np.asarray(rtmp), np.asarray(vtmp)

        # Vector difference and magnitudes
        dr, dv = r0 - rsgp4, v0 - vsgp4
        drm, dvm = np.linalg.norm(dr), np.linalg.norm(dv)
        
        # Exit check
        if (np.abs(drm) < 1e-1) and (np.abs(dvm) < 1e-3):
            converged = True
            break
        
        # Update state vector
        rnew = rnew + dr
        vnew = vnew + dv

    return newtle, converged

def classical_elements(r, v, epochyr, epochdoy, tle):
    """
    Classical elements from position and velocity. Following Chapter 4.4 from 
    'Methods of orbit determination for the micro computer' by Dan Boulet'.

    r: position vector in kilometers
    v: velocity vector in kilometers per second
    """

    # Constants
    xke = 0.07436680 # Gaussian gravitational constant
    xkmper = 6378.135 # Earth radius in kilometers
    ae = 1.0
    xmnpda = 1440.0 # Minutes per day
    mu = 1.0
    r2d = 180 / np.pi
    twopi = 2 * np.pi

    # Convert to fractional units
    r = np.asarray(r) / xkmper
    v = np.asarray(v) / (xke * xkmper * xmnpda / (ae * 86400))

    # Compute vector magnitudes and cross products
    rm = np.linalg.norm(r)
    vm2 = np.dot(v, v)
    rvm = np.dot(r, v)
    h = np.cross(r, v)
    chi = np.dot(h, h) / mu
    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    if (n[0] == 0) and (n[1] == 0):
        n[0] = 1
    nm = np.linalg.norm(n)
    
    # Eccentricity
    e = (vm2 / mu - 1 / rm) * r - rvm / mu * v
    ecc = np.linalg.norm(e)

    # Semi major axis
    a = (2 / rm - vm2 / mu)**(-1)

    # Inclination
    incl = np.arccos(h[2]/np.linalg.norm(h))

    # RA of ascending node
    node = np.mod(np.arctan2(n[1], n[0]), twopi)

    # Argument of perigee
    argp = np.mod(np.arccos(np.dot(n, e) / (nm * ecc)), twopi)
    if e[2] < 0:
        argp = twopi - argp

    # Mean anomaly
    xp = (chi - rm) / ecc
    yp = rvm / ecc * np.sqrt(chi / mu)
    b = a * np.sqrt(1 - ecc**2)
    cx = xp / a + ecc
    sx = yp / b
    ee = np.arctan2(sx, cx)
    m = np.mod(ee - ecc * sx, twopi)

    # Mean motion
    n = xke * np.sqrt(mu / a**3) * 1440 / twopi

    line0, line1, line2 = format_tle(tle.satno, epochyr, epochdoy, incl * r2d, node * r2d, ecc,
                                     argp * r2d, m * r2d, n, tle.bstar, tle.name, tle.desig, tle.classification,
                                     tle.ndot, tle.nddot, tle.ephtype, tle.elnum, tle.revnum)

    return TwoLineElement(line0, line1, line2)
