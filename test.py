#!/usr/bin/env python3
from twoline.twoline import TwoLineElement

if __name__ == "__main__":
    line0 = "0 ISS (ZARYA)"
    line1 = "1 25544U 98067A   20053.24093319  .00001847  00000-0  41427-4 0  9993"
    line2 = "2 25544  51.6429 202.1571 0004852 303.4083 121.5105 15.49190851214046"

    tle = TwoLineElement(line0, line1, line2)

    print(tle.name)
    print(tle.id)
    print(tle)
