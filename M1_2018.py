import BeamAnalysis
import numpy as np
def get_seconds_to_remove():
    """
    Returns seconds that are either extraneous data to the signal in question or are calibration seconds.
    """
    cal = []

    for x in range(1000):
        cal.append(x)
    for x in range(5150, 7201):
        cal.append(x)
    return cal

def frequencies():
    """
    A helper method that returns the median frequencies among the groups that were averaged together.
    """
    fre = []
    for i in range(512):
        fre.append(685.0 + i*0.244140625)
    return fre

timestamps = ['20180101234415','20180102004415', '20180102014415']


# add declination input
BeamAnalysis.main('C:/Python/BeamAnalysis','3srcNP', timestamps, frequencies(), 200,np.radians(22.01444),get_seconds_to_remove())
