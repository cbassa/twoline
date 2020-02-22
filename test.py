#!/usr/bin/env python3
from twoline.twoline import TwoLineElement

if __name__ == "__main__":
    line0 = "0 VANGUARD 3"
    line1 = "1 00020U 59007A   20052.88281644 -.00000025 +00000-0 -83416-5 0  9999"
    line2 = "2 00020 033.3418 237.3236 1666874 008.3442 354.1601 11.55682336217113"
    line0 = "IGS Radar 6"
    line1 = "1 43495U 18052A   19238.98700841 0.00000000  00000-0  00000-0 0    04"
    line2 = "2 43495  97.3707 353.7643 0001821 306.5562  53.4436 15.25981620    04"
    line0 = "0 ISS (ZARYA)"
    line1 = "1 25544U 98067A   20053.24093319  .00001847  00000-0  41427-4 0  9993"
    line2 = "2 25544  51.6429 202.1571 0004852 303.4083 121.5105 15.49190851214046"

    tle = TwoLineElement(line0, line1, line2)

    print(tle)

    print(line0)
    print(line1)
    print(line2)

    tle.format()
