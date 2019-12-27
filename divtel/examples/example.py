import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

print(np)

from divtel.telescope import Telescope, Array
from divtel import pointing

def hess_1():
    tel1 = Telescope(10*u.m, 0*u.m, 0*u.m, 20*u.m, 1*u.m)
    tel2 = Telescope(0 * u.m, 10 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    tel3 = Telescope(-10 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    tel4 = Telescope(0 * u.m, -10 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    return Array([tel1, tel2, tel3, tel4])


def hess_2():
    tel1 = Telescope(10*u.m, 0*u.m, 0*u.m, 20*u.m, 1*u.m)
    tel2 = Telescope(0 * u.m, 10 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    tel3 = Telescope(-10 * u.m, 0 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    tel4 = Telescope(0 * u.m, -10 * u.m, 0 * u.m, 20 * u.m, 1 * u.m)
    tel5 = Telescope(0 * u.m, 0 * u.m, 0 * u.m, 30 * u.m, 1 * u.m)
    return Array([tel1, tel2, tel3, tel4, tel5])


def random_array(n=10):
    tels = [Telescope(10*np.random.rand()*u.m,
                      10*np.random.rand() * u.m,
                      10*np.random.rand() * u.m,
                      np.random.rand() * u.m,
                      np.random.rand() * u.m,
                      )
            for i in range(n)
            ]
    return Array(tels)


def main():
    """
    just to test things
    """

    # print(tel1.position)
    # print(tel1.zenith)

    array = hess_1()
    # array = random_array()

    for tel in array.telescopes:
        tel.point_to_altaz(90*u.deg, 0*u.deg)

    # tel2.point_to_altaz(0*u.deg, -90*u.deg)
    # print(array.positions_array)
    # print(array.barycenter[0])
    #
    # ax = array.display_positions()
    # ax.legend()
    # plt.show()

    # print(array.pointing_vectors)

    alt = 40 * u.deg
    az = 0 * u.deg
    # alt = np.random.rand()*u.rad
    # az = np.random.rand()*u.rad

    array.divergent_pointing(0.1, alt, az)
    ax = array.display_2d(projection='xy')
    ax.legend()
    plt.show()

    print(array.positions_array)
    print(array.barycenter)


    ax = array.display_3d()
    ax.set_zlim(0, 1)
    plt.show()

    telescopes_distances = np.sqrt(np.sum((array.positions_array - array.barycenter)**2, axis=1))
    # print(telescopes_distances)
    tels_alt = np.array([tel.alt.value for tel in array.telescopes])
    tels_az = np.array([tel.az.value for tel in array.telescopes])
    # print(array.pointing_vectors)
    # print(tels_alt, tels_az)
    p = np.average(array.pointing_vectors, weights=telescopes_distances, axis=0)

    # print(np.average(tels_alt, weights=telescopes_distances), np.average(tels_az, weights=telescopes_distances))

    # print(p)
    # print(pointing.alt_az_to_vector(alt, az))
    # np.testing.assert_allclose(p, pointing.alt_az_to_vector(alt, az), rtol=1e-2)

    # from visualization import polar_stuff

    # tel1.point_to_altaz(70*u.deg, np.pi*u.rad)
    # ax = polar_stuff(plt.figure(), tel1)
    # ax.set_xlim(50, 80)
    # ax.set_ylim(50, 90)
    # plt.show()

    for tel in array.telescopes:
        print(tel.id)

if __name__ == "__main__":
    main()

