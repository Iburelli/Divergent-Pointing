import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt


class Telescope:

    def __init__(self, x, y, z, focal, camera_radius):

        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)
        self.alt = u.Quantity(0, u.rad)
        self.az = u.Quantity(0, u.rad)

    def point_to_altaz(self, alt, az):
        self.alt = alt.to(u.rad)
        self.az = az.to(u.rad)

    @property
    def zenith(self):
        return np.pi/2.*u.rad - self.alt

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

    @property
    def pointing_vector(self):
        return np.array([np.cos(self.alt)*np.cos(self.az), np.cos(self.alt)*np.sin(self.az), np.sin(self.az)])


class Array:

    def __init__(self, telescope_list):

        self.telescopes = {}
        for ii, tel in enumerate(telescope_list):
            self.telescopes[ii+1] = tel

    @property
    def positions_array(self):
        return np.array([tel.position for k, tel in self.telescopes.items()])

    @property
    def pointing_vectors(self):
        """
        all telescopes pointing vectors as an array

        Returns
        -------
        np.array
        """
        return np.array([tel.pointing_vector for key, tel in self.telescopes.items()])

    @property
    def barycenter(self):
        return self.positions_array.mean(axis=0)

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
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        ax.scatter(self.positions_array[:, 1], self.positions_array[:, 0], **kwargs, label='telescopes')
        ax.scatter(self.barycenter[0], self.barycenter[1], marker='+', label='barycenter')
        ax.quiver(self.positions_array[:, 1],
                  self.positions_array[:, 0],
                  self.pointing_vectors[:, 0],
                  self.pointing_vectors[:, 1],
                  color=kwargs['color']
                  )
        ax.set_ylabel('x [m]')
        ax.set_xlabel('y [m]')

        return ax

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
    for k, tel in array.telescopes.items():
        tel.point_to_altaz(70*u.deg, 0*u.deg)
    print(array.positions_array)
    print(array.barycenter[0])

    ax = array.display_positions()
    ax.legend()
    plt.show()

    print(array.pointing_vectors)

    # from visualization import polar_stuff

    # tel1.point_to_altaz(70*u.deg, np.pi*u.rad)
    # ax = polar_stuff(plt.figure(), tel1)
    # ax.set_xlim(50, 80)
    # ax.set_ylim(50, 90)
    # plt.show()


if __name__ == "__main__":
    main()