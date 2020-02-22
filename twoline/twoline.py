#!/usr/bin/env python3
import math
from datetime import datetime

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
            self.name = line0[2:]
        else:
            self.name = line0

        # Process line 1
        self.satno = int(line1[2:7])
        self.classification = line1[7]
        self.desig = line1[9:17].rstrip()
        self.epochyr = int(line1[18:20])
        self.epochdoy = float(line1[20:32])
        self.ndot = float(line1[33:43])
        fvalue = float(line1[44] + "." + line1[45:50])
        fexp = int(line1[50:52])
        self.nddot = fvalue * 10 ** fexp
        fvalue = float(line1[53] + "." + line1[54:59])
        fexp = int(line1[59:61])
        self.bstar = fvalue * 10 ** fexp
        self.ephtype = line1[62]
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
