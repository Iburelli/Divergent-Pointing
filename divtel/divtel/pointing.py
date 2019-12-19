"""
Functions to define telescopes pointings
We use the same reference frame as simtel_array:
X is pointing North
Y is pointing East
Z is pointing upward
Az is taken clock-wise from X (towards Y) and between -180 and 180 degrees
Alt is taken from ground (towards Z) and between -90 and 90 degrees
"""

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from telescope import Array, Telescope

def alt_az_to_vector(alt, az):
    """
    Compute a pointing vector from an alt,az pointing direction

    Parameters
    ----------
    alt: angle in rad
    az: angle in rad

    Returns
    -------
    np.array([x, y, z])
    """
    x = np.cos(alt) * np.cos(az)
    y = np.cos(alt) * np.sin(az)
    z = np.sin(alt)
    return np.array([x, y, z])


def retro_pointing(array, div, alt, az):

    B = array.barycenter
    norm = - np.log(div)
    Gx = B[0] - norm * np.cos(alt) * np.cos(az)
    Gy = B[1] - norm * np.cos(alt) * np.sin(az)
    Gz = B[2] - norm * np.sin(alt)
    return np.array([Gx, Gy, Gz])

def tel_div_pointing(tel, G):
    GT = np.sqrt(((tel.position - G) ** 2).sum())
    alt_tel = np.arcsin((tel.z.value - G[2]) / GT)
    az_tel = np.arctan2((tel.y.value - G[1]), (tel.x.value - G[0]))
    tel.point_to_altaz(alt_tel * u.rad, az_tel * u.rad)

def array_div_pointing(array, G):
    for ii, tel in array.telescopes.items():
        tel_div_pointing(tel, G)
        print(tel.alt, tel.az)

def main():
    tel1 = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    tel2 = Telescope(1 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    tel3 = Telescope(0 * u.m, 1 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)

    div = 1e-1
    alt = 70 * u.deg
    az = 90 * u.deg

    array = Array([tel1, tel2, tel3])
    for k, tel in array.telescopes.items():
        tel.point_to_altaz(alt, az)
    array.display_positions()
    plt.show()

    G = retro_pointing(array, div, alt, az)

    array_div_pointing(array, G)
    array.display_positions()
    plt.show()


if __name__ == "__main__":
    main()
