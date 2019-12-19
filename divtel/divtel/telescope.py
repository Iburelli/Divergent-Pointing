import numpy as np
import astropy.units as u
from astropy.coordinates import Angle


class Telescope:

    def __init__(self, x, y, z, focal, camera_radius):

        self.x = x.to(u.m)
        self.y = y.to(u.m)
        self.z = z.to(u.m)
        self.focal = focal.to(u.m)
        self.camera_radius = camera_radius.to(u.m)

    def point_to_altaz(self, alt, az):
        assert type(alt)==Angle
        assert type(az)==Angle
        self.alt = alt
        self.az = az


    def fov(self):
        """
        Area of the field of view in rad**2
        """
        return np.pi * (self.camera_radius / self.focal)**2



