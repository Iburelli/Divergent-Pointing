import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as plt


class Telescope:

    def __init__(self, x, y, z, focal, camera_radius):

        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = Angle(0*u.rad)
        self.az = Angle(0*u.rad)

    def point_to_altaz(self, alt, az):
        assert type(alt)==Angle
        assert type(az)==Angle
        self.alt = alt
        self.az = az

    @property
    def zenith(self):
        return Angle(np.pi/2.*u.rad) - self.alt

    @property
    def fov(self):
        """
        Area of the field of view in rad**2
        """
        return u.Quantity(np.pi * (self.camera_radius / self.focal)**2, u.rad**2)

    @property
    def position(self):
        return np.array([self.x.to(u.m).value, self.y.to(u.m).value, self.z.to(u.m).value]*u.m)

    def point_to_object(self):
        #TODO
        pass

class Array:

    def __init__(self, telescope_list):

        self.telescopes = {}
        for ii, tel in enumerate(telescope_list):
            self.telescopes[ii+1] = tel

    @property
    def positions_array(self):
        return np.array([tel.position for k, tel in self.telescopes.items()])

    def display_positions(self, ax=None, **kwargs):
        """
        Display the array

        Parameters
        ----------
        ax: `matplotlib.pyplot.axes`
        kwargs: args for `pyplot.scatter`

        Returns
        -------
        ax: `matplotlib.pyplot.axes`
        """
        ax = plt.gca() if ax is None else ax
        ax.scatter(self.positions_array[:,0], self.positions_array[:,1], **kwargs)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')

    def display_3d(self):
        #TODO: 3d representation of the array with the telescopes pointing
        pass


def main():
    """
    just to test things
    """
    tel1 = Telescope(0*u.m, 0*u.m, 0*u.m, 28*u.m, 1*u.m)
    tel2 = Telescope(1 * u.m, 0 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)
    tel3 = Telescope(0 * u.m, 1 * u.m, 0 * u.m, 28 * u.m, 1 * u.m)

    print(tel1.position)
    print(tel1.zenith)

    array = Array([tel1, tel2, tel3])
    print(array.positions_array)

    array.display_positions()
    plt.show()


if __name__ == "__main__":
    main()