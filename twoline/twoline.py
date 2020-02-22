#!/usr/bin/env python3

class TwoLineElement:
    """TLE class"""

    def __init__(self, line0, line1, line2):
        """Define a tle"""

        self.line0 = line0.rstrip()
        self.line1 = line1.rstrip()
        self.line2 = line2.rstrip()
        if self.line0[:2]=="0 ":
            self.name = self.line0[2:]
        else:
            self.name = self.line0
        self.id = self.line1.split(" ")[1][:5]

    def __repr__(self):
        return "%s\n%s\n%s"%(self.line0, self.line1, self.line2)
