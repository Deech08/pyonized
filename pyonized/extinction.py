from scipy.interpolate import interp1d

class Extinction:
    def __init__(self, dustmap, lbd):
        """
        Initialize the Extinction class with the given dustmap and coordinates

        Parameters
        ----------
        dustmap : object
            The dustmap query object ('Edenhofer2023Query', 'MarshallQuery', 'BayestarQuery').
        lbd : SkyCoord
            Sky coordinates for the extinction query.
        """
        
        self.dustmap = dustmap
        self.lbd = lbd
        self.extinction_curve = self._read_extinction_curve("extinction_curve.txt")
        self.extinction_data_ha, self.extinction_data_hb = self._get_extinction_parameters()

    def _read_extinction_curve(self, filepath):
        """
        Read the extinction curve from the file.

        Parameters
        ----------
        filepath : string
            Path to the file containing the extinction curve data.

        Returns
        -------
        dict : a dictionary with wavelength as keys (in nanometers) and extinction values as values
        """
        
        # Read the extinction curve from the file
        data = np.loadtxt(filepath, skiprows=1)
        wavelength = data[:, 0] # First column: wavelength in nanometers
        extinction = data[:, 1] # Second column: extinction values
        return {'wavelength': wavelength, 'extinction':extinction}
    
    def _get_extinction_parameters(self):
        """
        Get the extinction parameters for the specified dustmap.

        Returns
        -------
        extinction_data_ha, extinction_data_hb : array
            Extinction data for H-alpha and H-beta.
        """

        if self.dustmap is None:
            return np.ones_like(self.lbd.l.value), np.ones_like(self.lbd.l.value) # Default extinction values if dustmap is None

        if isinstance(self.dustmap, Edenhofer2023Query):
            eden_extinction = self.dustmap(self.lbd, mode='mean')
            
            extinction_data_ha = self._convert_zgr23_to_extinction(656.3) * eden_extinction * 2.8
            extinction_data_hb = self._convert_zgr23_to_extinction(486.1) * eden_extinction * 2.8

        elif isinstance(self.dustmap, MarshallQuery):
            wave_Ks = 2.17 * u.micron
            A_KS_to_A_v = 1. / extinction_law(np.array([wave_Ks.to(u.AA).value]), 1.)
            wave_ha = np.array([6562.8]) * u.AA
            A_V_to_A_ha = extinction_law(wave_ha.to(u.AA).value, 1.)
            wave_hb = np.array([4861.32]) * u.AA
            A_V_to_A_hb = extinction_law(wave_hb.to(u.AA).value, 1.)
            
            extinction_data_ha = 10 ** (-0.4 * A_KS_to_A_v * A_V_to_A_ha * self.dustmap.query(self.lbd))
            extinction_data_hb = 10 ** (-0.4 * A_KS_to_A_v * A_V_to_A_hb * self.dustmap.query(self.lbd))

        elif isinstance(self.dustmap, BayestarQuery):
            bayestar_extinction = self.dustmap.query(self.lbd)
            
            extinction_data_ha = self._convert_zgr23_to_extinction(656.3) * bayestar_extinction * 2.742
            extinction_data_hb = self._convert_zgr23_to_extinction(486.1) * bayestar_extinction * 2.742

        else:
            raise ValueError(f"Unknown dustmap class: {dustmap}")

        return extinction_data_ha, extinction_data_hb

    def _convert_zgr23_to_extinction(self, wavelength):
        """
        Convert the ZGR23 extinction to the specified wavelength.

        Parameters
        ----------
        zgr23_extinction : float
            ZGR23 extinction value.
        wavelength : float
            Wavelength in nanometers.

        Returns
        -------
        float
            Extinction value for the specified wavelength.
        """
        
        wavelengths = self.extinction_curve['wavelength']
        extinctions = self.extinction_curve['extinction']
        interp_func = interp1d(wavelengths, extinctions, kind='linear', fill_value="extrapolate")
        extinct = interp_func(wavelength)
        return extinct
        