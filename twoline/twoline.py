#!/usr/bin/env python3
import math
import numpy as np
from datetime import datetime
from astropy.time import Time
from sgp4.api import Satrec, jday

# Using definition from https://celestrak.com/columns/v04n03/
class TwoLineElement:
    """TLE class"""

    def __init__(self, line0, line1, line2):
        """Define a tle"""

        # Clean input lines
        self.line0 = line0.rstrip()
        self.line1 = line1.rstrip()
        self.line2 = line2.rstrip()

        # Extract object name
        if self.line0[:2] == "0 ":
            self.name = self.line0[2:]
        else:
            self.name = self.line0

        # Process line 1
        self.satno = int(line1[2:7])
        self.classification = line1[7]
        self.desig = line1[9:17].rstrip()
        self.desig_year = line1[9:11]
        self.desig_id = line1[11:17].rstrip()
        self.epochyr = int(line1[18:20])
        self.epochdoy = float(line1[20:32])
        self.ndot = float(line1[33:43])
        fvalue = float(line1[44] + "." + line1[45:50])
        fexp = int(line1[50:52])
        self.nddot = fvalue * 10 ** fexp
        fvalue = float(line1[53] + "." + line1[54:59])
        fexp = int(line1[59:61])
        self.bstar = fvalue * 10 ** fexp
        self.ephtype = int(line1[62])
        self.elnum = int(line1[64:68])
        
        # Process line 2
        satno = int(line2[2:7])
        self.incl = float(line2[8:16])
        self.node = float(line2[17:25])
        self.ecc = float("0." + line2[26:33].replace(" ", "0"))
        self.argp = float(line2[34:42])
        self.m = float(line2[43:51])
        self.n = float(line2[52:63])
        self.revnum = int(line2[63:68])

        # Compute epoch
        self.epoch = epoch_to_datetime(self.epochyr, self.epochdoy)
   
    def __repr__(self):
        return "%s\n%s\n%s"%(self.line0, self.line1, self.line2)

    def __str__(self):
        text = f"Satellite name:                        {self.name:>12}\n" \
               f"Satellite number:                      {self.satno:>12}\n" \
               f"International designator:              {self.desig:>12}\n" \
               f"Classification:                        {self.classification:>12}\n" \
               f"Epoch:                               {self.epochyr:02d}{self.epochdoy:012.8f}\n" \
               f"Epoch (datetime):        {self.epoch}\n" \
               f"First time derivative of mean motion:  {self.ndot:>12}\n" \
               f"Second time derivative of mean motion: {self.nddot:>12}\n" \
               f"BSTAR drag term:                       {self.bstar:>12}\n" \
               f"Ephemeris type:                        {self.ephtype:>12}\n" \
               f"Element number:                        {self.elnum:>12}\n" \
               f"Inclination (deg):                     {self.incl:>12}\n" \
               f"RA of the ascending node (deg):        {self.node:>12}\n" \
               f"Eccentricity:                          {self.ecc:>12}\n" \
               f"Argument of perigee (deg):             {self.argp:>12}\n" \
               f"Mean anomaly (deg):                    {self.m:>12}\n" \
               f"Mean motion (revolutions/day):         {self.n:>12}\n" \
               f"Revolution number at epoch:            {self.revnum:>12}"         
        return text

    def print_tle(self):
        print(f"{self.line0}\n{self.line1}\n{self.line2}")
    
    def format(self):
        return format_tle(self.satno, self.epochyr, self.epochdoy, self.incl, self.node, self.ecc, self.argp, self.m, self.n,
                          self.bstar, self.name, self.desig, self.classification,
                          self.ndot, self.nddot, self.ephtype, self.elnum, self.revnum)

    def from_parameters(satno, epochyr, epochdoy, incl, node, ecc, argp, m, n, bstar=0, name="Test object", desig="20500A", classification="U",
                        ndot=0, nddot=0, ephtype=0, elnum=0, revnum=0):

        
        # Compute epoch
        epoch = epoch_to_datetime(epochyr, epochdoy)

        # Format ndot term
        ndotstr = format_ndot(ndot)
        
        # Format nddot term
        nddotstr = format_bstar_or_nddot(nddot)
        
        # Format B* drag term
        bstarstr = format_bstar_or_nddot(bstar)
        
        # Format lines
        line0 = f"{name}"
        line1 = set_checksum(f"1 {satno:05d}{classification:1s} {desig:8s} {epochyr:02d}{epochdoy:012.8f} " \
            f"{ndotstr:10s} {nddotstr:8s} {bstarstr:8s} {ephtype:1d} {elnum:4d}0")
        line2 = set_checksum(f"2 {satno:05d} {incl:8.4f} {node:8.4f} {ecc * 10000000:07.0f} {argp:8.4f} " \
            f"{m:8.4f} {n:11.8f}{revnum:5d}0")

        return TwoLineElement(line0, line1, line2)
    
    def propagate(self, tnew, drmin=1e-1, dvmin=1e-3, niter=100):
        # New epoch
        epochyr, epochdoy = datetime_to_epoch(tnew)
        jd, fr = jday(tnew.year, tnew.month, tnew.day, tnew.hour, tnew.minute, tnew.second + tnew.microsecond / 1000000)
    
        # Compute reference state vector
        sat = Satrec.twoline2rv(self.line1, self.line2)
        e, r, v = sat.sgp4(jd, fr)
        r0, v0 = np.asarray(r), np.asarray(v)

        # Start loop
        rnew, vnew = r0, v0
        converged = False
        for i in range(niter):
            # Convert state vector into new TLE
            newtle = classical_elements(rnew, vnew, epochyr, epochdoy, self)
        
            # Propagate with SGP4
            sat = Satrec.twoline2rv(newtle.line1, newtle.line2)
            e, rtmp, vtmp = sat.sgp4(sat.jdsatepoch, sat.jdsatepochF)
            rsgp4, vsgp4 = np.asarray(rtmp), np.asarray(vtmp)

            # Vector difference and magnitudes
            dr, dv = r0 - rsgp4, v0 - vsgp4
            drm, dvm = np.linalg.norm(dr), np.linalg.norm(dv)
        
            # Exit check
            if (np.abs(drm) < drmin) and (np.abs(dvm) < dvmin):
                converged = True
                break
        
            # Update state vector
            rnew = rnew + dr
            vnew = vnew + dv

        return newtle, converged


# Format ndot term
def format_ndot(x):
    if x < 0:
        sign = "-"
    else:
        sign = " "
    value = f"{math.fabs(x):10.8f}"
    return f"{sign:1s}{value[1:]:9s}"
    
# Format bstar or nddot drag
def format_bstar_or_nddot(x):
    if x == 0.0:
        return " 00000-0"
    else:
        fexp = int(math.log(math.fabs(x)) / math.log(10))
        fvalue = (math.fabs(x) / 10 ** fexp) * 100000
        if x < 0:
            sign = "-"
        else:
            sign = " "
        return f"{sign:1s}{fvalue:05.0f}{fexp:+02d}"

    
# Compute checksum
def set_checksum(line):
    s = 0
    for c in line[:-1]:
        if c.isdigit():
            s += int(c)
        if c == "-":
            s += 1
    return line[:-1]+"%d"%(s%10)        
    
# From Jean Meeus' Astronomical Algorithms (chapter 7)
def epoch_to_datetime(epochyr, epochdoy):
    # Determine year (will fail after 2056!)
    if epochyr < 57:
        year = epochyr + 2000
    else:
        year = epochyr + 1900

    # Leap year?
    if (year % 4 == 0) and (year % 400 !=0):
        k = 1
    else:
        k = 2

    # Find month
    n = math.floor(epochdoy)
    month = math.floor(9 * (k + n) / 275 + 0.98)
    if n < 32:
        month = 1

    # Find day
    day = n - math.floor(275 * month / 9) + k * math.floor((month + 9) / 12) + 30

    # Find time
    x = 24 * (epochdoy - n)
    hour = math.floor(x)
    x = 60 * (x - hour)
    minute = math.floor(x)
    second = 60 * (x - minute)
    isecond = math.floor(second)
    usecond = math.floor(1000000 * (second - isecond))

    return datetime(year, month, day, hour, minute, isecond, usecond)

# From Jean Meeus' Astronomical Algorithms (chapter 7)
def datetime_to_epoch(t):
    # Extract information
    year, month, day, hour, minute, isecond, usecond = t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond

    # Leap year?
    if (year % 4 == 0) and (year % 400 !=0):
        k = 1
    else:
        k = 2

    # Day of year
    n = math.floor(275 * month / 9) - k * math.floor((month + 9) / 12) + day - 30

    # Fractional day
    fday = hour / 24 + minute / 1440 + (isecond + usecond / 1000000) / 86400
    
    # Determine year (will fail after 2056!)
    if year < 2000:
        epochyr = year - 1900
    else:
        epochyr = year - 2000

    epochdoy = n + fday

    return epochyr, epochdoy


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

