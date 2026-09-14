# This code modified based on repo https://github.com/jrleeman/Crust1.0.git

import os
from math import floor

import numpy as np


class CrustModel:
    """Load the nine-layer CRUST1.0 global one-degree model.

    Parameters
    ----------
    path_crust1 : str
        Directory containing crust1.vp, crust1.vs, crust1.rho and crust1.bnds.

    Returns
    -------
    model : CrustModel
        Model arrays vp, vs, rho and bnds with shape (180, 360, 9).

    Raises
    ------
    OSError
        One of the four CRUST1 files cannot be read.
    ValueError
        The model files do not contain the expected grid size.

    Notes
    -----
    Vp/Vs are km/s, density is g/cm3, and layer boundary elevations are km relative to sea level. Only get_point is part of the supported query interface.
    """

    def __init__(self, path_crust1):
        # Read in data files
        self.vp = np.loadtxt(os.path.join(path_crust1, "crust1.vp"))
        self.vs = np.loadtxt(os.path.join(path_crust1, "crust1.vs"))
        self.rho = np.loadtxt(os.path.join(path_crust1, "crust1.rho"))
        self.bnds = np.loadtxt(os.path.join(path_crust1, "crust1.bnds"))

        # Reshape to a lon,lat,layer grid. The 0,0 index value
        # is at 90 south and 180 latitude.
        self.vp = self.vp.reshape((180, 360, 9))
        self.vs = self.vs.reshape((180, 360, 9))
        self.rho = self.rho.reshape((180, 360, 9))
        self.bnds = self.bnds.reshape((180, 360, 9))

        self.layer_names = [
            "water",
            "ice",
            "upper_sediments",
            "middle_sediments",
            "lower_sediments",
            "upper_crust",
            "middle_crust",
            "lower_crust",
            "mantle",
        ]

    def _get_index(self, lat, lon):
        """
        Returns in index values used to query the model for a given lat lon.

        Paramaters
        ----------
        lat : float
        Latitude of interest

        lat : flaot
        Longitude of interest

        Returns
        -------
        ilat : int
        Index for given latitude

        ilon : int
        Index for given longitude
        """

        # Make sure the longitude is between -180 and 180
        if lon > 180:
            lon -= 360
        if lon < -180:
            lon += 360

        # Find the index in the data for given lat and lon.
        # The cells are 1 deg wide, so the poles and lon = 180 land one past the
        # last cell and have to be clamped back onto it.
        ilat = min(179, max(0, floor(90.0 - lat)))
        ilon = min(359, max(0, floor(180 + lon)))

        return int(ilat), int(ilon)

    def get_point(self, lat, lon):
        """Select the CRUST1.0 grid cell at a geographic location.

        Parameters
        ----------
        lat : float
            Latitude in degrees.
        lon : float
            Longitude in degrees.

        Returns
        -------
        layers : dict
            Layer names map to [Vp km/s, Vs km/s, density g/cm3, thickness km,
            top elevation km]. Layers thinner than 0.01 km are omitted except mantle.

        Notes
        -----
        Grid selection is nearest enclosing one-degree cell, without spatial interpolation. Depth is the negative of elevation.
        """

        # Get index for arrays of data at this location
        ilat, ilon = self._get_index(lat, lon)

        # Calculate the thickness of the layers, add zero to the end
        # for the mantle since it's not defined
        thickness = np.abs(np.ediff1d(self.bnds[ilat, ilon], to_end=[0]))

        model_layers = dict()

        for i, layer in enumerate(self.layer_names):
            vp = self.vp[ilat, ilon][i]
            vs = self.vs[ilat, ilon][i]
            rho = self.rho[ilat, ilon][i]
            bnd = self.bnds[ilat, ilon][i]
            layer_thickness = thickness[i]

            # If the layer has thickness or is the mantle, write it
            if layer_thickness >= 0.01 or layer == "mantle":
                model_layers[layer] = [vp, vs, rho, layer_thickness, bnd]

        return model_layers
