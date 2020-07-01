#!/usr/bin/env python3
from twoline.twoline import TwoLineElement
from datetime import datetime

if __name__ == "__main__":
    # Input TLE
    line0 = "0 VANGUARD 3"
    line1 = "1 00020U 59007A   20052.88281644 -.00000025 +00000-0 -83416-5 0  9999"
    line2 = "2 00020 033.3418 237.3236 1666874 008.3442 354.1601 11.55682336217113"

    # Read TLE
    tle = TwoLineElement(line0, line1, line2)
    tle.print_tle()

    # Set new date
    tnew = datetime(2020, 2, 23, 1, 2, 4, 0)

    # Propagate TLE to new date
    newtle, converged = tle.propagate(tnew)
    newtle.print_tle()

    # Print propagation state
    print(converged)

    # Define manual TLE
    tle = TwoLineElement.from_parameters(99999, 20, 086.12345, 97.2, 100.0, 0.01234, 123.4567, 0123.4567, 14.567890, 5e-5)
    tle.print_tle()
