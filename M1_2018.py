import BeamAnalysis

def get_seconds_to_remove():
    """
    Returns seconds that are either extraneous data to the signal in question or are calibration seconds.
    """
    cal = []
    for i in range(30):
        for j in range(3):
            cal.append(208+j+240*i)
    # these are data that are extraneous to M1 signal

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

def baseline_to_feeds():
    """
    Returns two lists with corresponding Chnno -> feed no. pairs.
    """
    Chnfeed = open('C:/Python/BeamAnalysis/dish_bl_order.txt', 'r')
    lines=list(Chnfeed)
    feeds = []
    Chnnos = []
    for val in range(0, len(lines), 2):
        lines[val]=str(lines[val][1:-3])
        lines[val] = lines[val].split()
        for i in range(len(lines[val])):
            lines[val][i]=int(lines[val][i])
        feeds.append(lines[val])
    for val in range(1, len(lines), 2):
        lines[val]=int(str(lines[val][:-2]))
        Chnnos.append(lines[val])
    return feeds,Chnnos
feeds, chnnos = baseline_to_feeds()
chntofe=[feeds,chnnos]

BeamAnalysis.main('C:/Python/BeamAnalysis','3srcNP', timestamps, frequencies(), chntofe, get_seconds_to_remove())

