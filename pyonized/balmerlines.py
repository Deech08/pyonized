import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord

import dask.array as da
import pyneb as pn
from scipy import integrate
from scipy.integrate import simpson
from scipy.interpolate import interp1d

from extinction import fm07 as extinction_law
from dustmaps.config import config
from dustmaps.edenhofer2023 import Edenhofer2023Query
from dustmaps.marshall import MarshallQuery
from dustmaps.bayestar import BayestarQuery

class BalmerLines:
    def __init__(self, ncells, maxdist, dustmap_class=None, l=6*u.deg, b=5*u.deg, dl=None, sigma_i=None, T_e=None, velocities=None): 
        """
        This class defines length of arrays, distance between cells, distances from zero,
        electron velocity, electron density, sky coordinates, and accounts for dust by
        calculating the extinction.

        Parameters
        ----------
        self : the instance of the class
        ncells : number of cells/inputs in the array
        maxdist : maximum distance
        dustmap_class : class, optional
            The dustmap query class to use (Edenhofer2023Query, MarshallQuery, or BayestarQuery). Default is None.
        l : float, optional
            Galactic longitude in degrees. Default is 6.
        b : float, optional
            Galactic latitude in degrees. Default is 5.
        dl : float with units, optional
            Distance between cells. Default is None.
        sigma_i : float with units, optional
            Electron velocity dispersion. Default is None.
        T_e : float with units, optional
            Electron temperature. Default is None.
        velocities : float with units, optional
            Electron velocity. Default is None.
        """
        
        self.ncells = ncells
        self.maxdist = maxdist * u.kpc
        self.dl = maxdist/ncells * u.pc
        self.distance = np.arange(0, maxdist.to_value(), maxdist.to_value() / ncells) * u.kpc
        self.v_e = np.zeros(ncells) * (u.km / u.second)
        self.n_e = np.zeros(ncells) / u.cm**(-3)
        self.l = np.ones(ncells) * l
        self.b = np.ones(ncells) * b
        self.dustmap_class = dustmap_class

        # Sky Coordinates (Latitude, Longitude, Distance)
        lbd = SkyCoord(l = self.l ,
                       b = self.b,
                       distance = self.distance,
                       frame = "galactic")

        # Get extinction data and parameters based on the specified dustmap
        if dustmap_class is BayestarQuery:
            dustmap = dustmap_class(max_samples=1)
        elif dustmap_class is not None:
            dustmap = dustmap_class()
        else:
            dustmap = None
        
        extinction = Extinction(dustmap, lbd)
            
        # H-alpha
        self.trans_ha = extinction.extinction_data_ha
        
        # H-beta
        self.trans_hb = extinction.extinction_data_hb
        
        # Default values
        self.dl = dl if dl is not None else maxdist/ncells
        self.sigma_i = sigma_i if sigma_i is not None else 12 * (u.km / u.second)
        self.T_e = T_e if T_e is not None else 10000 * u.K
        self.velocities = velocities if velocities is not None else np.linspace(-300, 300, 1000)
        
        # Precompute intensities
        self.precompute_intensities()
    
    def manual_input(self, index, v_e_value, n_e_value):
        """
        This function allows for manual input into the lists of electron velocity and density.

        Parameters
        ----------
        self : the instance of the class
        index : int
            The index at which the new value will be inserted.
        v_e_value : float
            The new electron velocity value to be assigned at the specified index.
        n_e_value : float
            The new electron density value to be assigned at the specified index.
        """
        
        self.v_e[index] = v_e_value * (u.km / u.second)
        self.n_e[index] = n_e_value / u.cm**(-3)
        self.precompute_intensities()
    
    def _calculate_intensity(self, velocities, const, trans, b_lambda, extinction_correction):
        """
        This function calculates the intensity of emission lines at a range of velocities.
        The assumptions for each variable in this function is based on Krishnarao, Benjamin,
        and Haffner (2020).

        Parameters
        ----------
        self : the instance of the class
        velocities : array
            Array of input velocities
        const : intensity constant of an emission line in units: R/kms^-1
        trans : array
            Extinction-corrected transmission value for each cell
        b_lambda : constant for a specific emission line based on Draine's "Physics of the Instellar and Intergalatic Medium" 

        Returns
        -------
        intensities : array
            Array of intensities for emission lines.
        """
        
        velocities = velocities * (u.km / u.second)
        dimensionless_T = (self.T_e / (10**4 * u.K))
        b_lambda_factor = b_lambda
        dimensionless_vel = (-0.5 * ((velocities[:, np.newaxis] - self.v_e) / self.sigma_i) ** 2).to_value()
        I = const * (self.n_e**2) * self.distance * (self.sigma_i**-1) * (dimensionless_T ** b_lambda_factor) * np.exp(dimensionless_vel)
        if extinction_correction:
            I = I * trans
        return np.sum(I, axis=1)

    def calculate_intensity_ha(self, velocities, extinction_correction):
        """
        This function calculates the intensity of H-alpha emission line at a range of velocities.

        Parameters
        ----------
        self : the instance of the class
        velocities : array
            Array of input velocities.

        Returns
        -------
        intensities : array
            Array of intensities for H-alpha emission line.
        """
        
        H1 = pn.RecAtom('H', 1)
        emissivity_ha = (H1.getEmissivity(self.T_e.value, self.n_e.value, lev_i=3, lev_j=2))
        const_ha = ((emissivity_ha*(3.09*(10**18)))/((4 * np.pi * u.steradian) * np.sqrt(2 * np.pi))) * (1 / (3.026*(10**-12))) * (1 / ((10**6) / (4 * np.pi))) 
        b_lambda_ha = -0.942 - 0.031 * np.log(dimensionless_T)
        return self._calculate_intensity(velocities, const_ha, self.trans_ha, b_lambda_ha, extinction_correction)

    def calculate_intensity_hb(self, velocities, extinction_correction):
        """
        This function calculates the intensity of H-beta emission line at a range of velocities.

        Parameters
        ----------
        self : the instance of the class
        velocities : array
            Array of input velocities.

        Returns
        -------
        intensities : array
            Array of intensities for H-beta emission line.
        """
        
        H1 = pn.RecAtom('H', 1)
        emissivity_hb = (H1.getEmissivity(self.T_e.value, self.n_e.value, lev_i=4, lev_j=2))
        const_hb = ((emissivity_hb*(3.09*(10**18)))/((4 * np.pi * u.steradian) * np.sqrt(2 * np.pi))) * (1 / (3.026*(10**-12))) * (1 / ((10**6) / (4 * np.pi)))
        b_lambda_hb = -0.874 - 0.058 * np.log(dimensionless_T)
        return self._calculate_intensity(velocities, const_hb, self.trans_hb, b_lambda_hb, extinction_correction)

    def precompute_intensities(self):
        """
        This function precomputes the intensities of emission lines at a range
        of velocities, taking into account whether there is extinction correction.

        Parameters
        ----------
        self : the instance of the class
        extinction_correction : boolean
            Flag indicating whether extinction correction is applied or not.
        """
        
        self.intensities_ha_with = self.calculate_intensity_ha(self.velocities, extinction_correction=True)
        self.intensities_hb_with = self.calculate_intensity_hb(self.velocities, extinction_correction=True)
        self.intensities_ha_without = self.calculate_intensity_ha(self.velocities, extinction_correction=False)
        self.intensities_hb_without = self.calculate_intensity_hb(self.velocities, extinction_correction=False)

    def plot_intensity(self, line='ha', extinction_correction=True):
        """
        This function plots individual graphs of the intensities for H-alpha and H-beta
        emission lines, both with and without extinction correction.

        Parameters
        ----------
        self : the instance of the class
        line : indicates the type of emission line
        extinction_correction : boolean
            Flag indicating whether extinction correction is applied or not. Default is True.
        """
        
        if line == 'ha':
            intensities = self.intensities_ha_with if extinction_correction else self.intensities_ha_without
            label = 'H-alpha Intensity'
        else:
            intensities = self.intensities_hb_with if extinction_correction else self.intensities_hb_without
            label = 'H-beta Intensity'
            color = 'orange'
        
        plt.plot(self.velocities, intensities, label=label)
        plt.title(f'{label} {"with" if extinction_correction else "without"} dust correction for multiple cells at a range of velocities')
        plt.xlabel('Input velocity (km/s)')
        plt.ylabel(f'I_{label} (R/(km s^-1))')
        plt.legend()
        plt.show()

    def plot_comparison(self, velocities):
        """
        This function plots the comparison of intensities for H-alpha and H-beta emission lines,
        both with and without extinction correction.

        Parameters
        ----------
        self : the instance of the class
        velocities : array
            Array of input velocities.
        """
        
        plt.plot(self.velocities, self.intensities_ha_with, label='H-alpha Intensity with Dust Correction', color='blue')
        plt.plot(self.velocities, self.intensities_hb_with, label='H-beta Intensity with Dust Correction', color='orange')
        plt.plot(self.velocities, self.intensities_ha_without, label='H-alpha Intensity without Dust Correction', linestyle='dashed', color='blue')
        plt.plot(self.velocities, self.intensities_hb_without, label='H-beta Intensity without Dust Correction', linestyle='dashed', color='orange')
        plt.title('Comparison of H-alpha and H-beta Intensities with and without dust correction')
        plt.xlabel('Input velocity (km/s)')
        plt.ylabel('Intensity (R/(km s^-1))')
        plt.legend()
        plt.show()

    def plot_integrated_spectra(self, extinction_correction=True):
        """
        This function plots the integrated spectra for H-alpha and H-beta emission lines,
        both with and without extinction correction.

        Parameters
        ----------
        self : the instance of the class
        extinction_correction : boolean
            Flag indicating whether extinction correction is applied or not. Default is True.
        """
        
        wavelengths = np.linspace(400, 700, 1000)  # Range of wavelengths from 400 nm to 700 nm

        # Gaussian profile function
        def gaussian(x, amp, mean, stddev):
            return amp * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

        # Calculate integrated intensities
        if extinction_correction:
            ha_peak_intensity = np.max(self.intensities_ha_with.value) # .value strips units
            hb_peak_intensity = np.max(self.intensities_hb_with.value)
            label = 'with Dust Correction'
        else:
            ha_peak_intensity = np.max(self.intensities_ha_without.value)
            hb_peak_intensity = np.max(self.intensities_hb_without.value)
            label = 'without Dust Correction'

        # Create Gaussian profiles for H-alpha and H-beta
        h_alpha_params = {'amp': ha_peak_intensity, 'mean': 656.3, 'stddev': 1.0}
        h_beta_params = {'amp': hb_peak_intensity, 'mean': 486.1, 'stddev': 1.0}

        spectra_h_alpha = gaussian(wavelengths, **h_alpha_params)
        spectra_h_beta = gaussian(wavelengths, **h_beta_params)
        total_spectra = spectra_h_alpha + spectra_h_beta

        # Integrating the spectral lines for display purposes
        integral_h_alpha = simpson(spectra_h_alpha, wavelengths)
        integral_h_beta = simpson(spectra_h_beta, wavelengths)

        # Plotting the results
        plt.figure(figsize=(10, 6))
        plt.plot(wavelengths, total_spectra, label='Total Spectra')
        plt.plot(wavelengths, spectra_h_alpha, label='H-alpha')
        plt.plot(wavelengths, spectra_h_beta, label='H-beta')
        plt.fill_between(wavelengths, spectra_h_alpha, alpha=0.3)
        plt.fill_between(wavelengths, spectra_h_beta, alpha=0.3)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity (R)')
        plt.title(f'Integrated Spectra of H-alpha and H-beta {label}')
        plt.legend()
        plt.grid(True)

        # Showing the integral results
        plt.text(450, 0.5, f'H-alpha Integral: {integral_h_alpha:.2f}', fontsize=12, color='blue')
        plt.text(450, 0.4, f'H-beta Integral: {integral_h_beta:.2f}', fontsize=12, color='orange')

        plt.show()


