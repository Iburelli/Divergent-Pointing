import numpy as np
import astropy.units as u
from astropy.coordinates import Angle

from telescope import Array, Telescope

def alt_az_to_vector(alt, az):
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
    GT = np.sqrt(((tel - G) ** 2).sum())
    alt_tel = np.arcsin((tel[2] - G[2]) / GT)
    az_tel = np.arctan2((tel[1] - G[1]), (tel[0] - G[0]))
    #         print("alt", np.rad2deg(alt_tel))
    #         print("az", az_tel)
    return 90 - np.rad2deg(alt_tel), np.rad2deg(np.mod(az_tel, 2 * np.pi))

def array_div_pointing(array, G):
    pointings = []
    for ii, tel in enumerate(array.positions_array):
        theta, phi = tel_div_pointing(tel, G)
        print("theta: ", theta)
        print("phi: ", phi)

def main():
    tel1 = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    tel2 = Telescope(1 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    tel3 = Telescope(0 * u.m, 1 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)

    div = 1e-15
    alt = np.deg2rad(90)
    az = np.pi

    array = Array([tel1, tel2, tel3])
    G = retro_pointing(array, div, alt, az)

    array_div_pointing(array, G)
    

if __name__ == "__main__":
    main()
