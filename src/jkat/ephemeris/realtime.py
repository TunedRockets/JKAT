''' 
file for (shudders in horror) calculating time.
this is complicated and horrible, but needs to be done

"natural" time in use is with the epoch J2000 in TT
that is, time zero is
Gregorian date jan 1, 2000 at 12:00.
which is equivalent to:
- julian date 2451545.0 TT
- jan 1, 2000, 11:59:27.816 TAI
- jan 1, 2000, 11:58:55.816 UTC
- 0.0 MJD2000
- 2_451_545.0 JD
- something else in other JD systems

assume datetime is UTC
'''
import datetime as dt
from ..utils.consts import DAY

__all__ = [
    'to_time',
    'from_time',
    'J2000_epoch',
    'to_JD',
    'from_JD'
]

J2000_epoch = dt.datetime(
    2000, 1, 1,
    11, 59, 27, 816_000
)
'''epoch of the J2000 / ICRF timekeeping system (in UTC)'''

# todo for now:
# turn python datetime into float time from epoch
# reverse this
# turn simple date into MJD-TT
# in future:
# be more intelligent with different conversion


def to_time(d:dt.datetime|dt.date)->float:
    '''turn python datetime (in UTC) into seconds since J2000 epoch'''
    if isinstance(d,dt.date): d = dt.datetime(d.year,d.month,d.day)
    delta = d - J2000_epoch
    return delta.total_seconds()


def from_time(t:float)->dt.datetime:
    delta = dt.timedelta(seconds=t)
    return J2000_epoch + delta

def to_JD(t:float, MJD2000:bool=False)->float:
    '''return julian date from time since ephemeris,
    if MJD2000 is true will use 2000 as 0 JD'''

    t /= DAY
    if not MJD2000: t += 2451545
    return t

def from_JD(JD:float, MJD2000:bool=False)->float:
    '''return time since ephemeris from julian date,
    if MJD2000 is true will use 2000 as 0 JD'''

    if not MJD2000: JD -= 2451545
    JD *= DAY
    return JD
