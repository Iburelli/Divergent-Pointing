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
        return np.array([np.cos(self.alt.to(u.rad))*np.cos(self.az.to(u.rad)),
                         np.cos(self.alt.to(u.rad))*np.sin(self.az.to(u.rad)),
                         np.sin(self.alt.to(u.rad))])


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

    def display_2d(self, projection='xy', ax=None, **kwargs):
        """
        Display the array

        Parameters
        ----------
        projection: str
            'xy', 'xz' or 'yz'
        ax: `matplotlib.pyplot.axes`
        kwargs: args for `pyplot.scatter`

        Returns
        -------
        ax: `matplotlib.pyplot.axes`
        """
        ax = plt.gca() if ax is None else ax
        if 'color' not in kwargs:
            kwargs['color'] = 'black'

        if projection=='xy':
            xx = self.positions_array[:, 1]
            yy = self.positions_array[:, 0]
            xb = self.barycenter[1]
            yb = self.barycenter[0]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 0]
            xlabel = 'y [m]'
            ylabel = 'x [m]'

        elif projection=='xz':
            xx = self.positions_array[:, 0]
            yy = self.positions_array[:, 2]
            xb = self.barycenter[0]
            yb = self.barycenter[2]
            xv = self.pointing_vectors[:, 0]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'x [m]'
            ylabel = 'z [m]'

        elif projection=='yz':
            xx = self.positions_array[:, 1]
            yy = self.positions_array[:, 2]
            xb = self.barycenter[1]
            yb = self.barycenter[2]
            xv = self.pointing_vectors[:, 1]
            yv = self.pointing_vectors[:, 2]
            xlabel = 'y [m]'
            ylabel = 'z [m]'

        else:
            breakpoint()


        ax.scatter(xx, yy, **kwargs, label='telescopes')
        ax.scatter(xb, yb, marker='+', label='barycenter')
        ax.quiver(xx, yy, xv, yv,
                  color=kwargs['color']
                  )
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.grid('on')
        ax.axis('equal')

        return ax

    def display_3d(self):
        #TODO: 3d representation of the array with the telescopes pointing
        pass



