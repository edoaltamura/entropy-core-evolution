.. _profiles:

Radial profiles
===============

We compute three radial profiles: density, mass-weighted temperature, and entropy. They are
centred in the centre of potential of the main object, and are computed from the particle
(-weighted) averages in 30 log-spaced spherical shells between :math:`0.05~r_{500}` to
:math:`2.5~r_{500}`.


The logarithmic centre of the radial bins is provided and can be loaded as

>>> radial_bin_centers = np.load("data/VR2915_+1res_ref/radial_bin_centers.npy")
>>> radial_bin_centers.shape
(30,)
>>> radial_bin_centers
array([0.00554565, 0.00682211, 0.00839236, 0.01032404, 0.01270034,
       0.0156236 , 0.01921971, 0.02364354, 0.02908561, 0.03578029,
       0.04401589, 0.0541471 , 0.06661022, 0.08194199, 0.1008027 ,
       0.12400461, 0.15254694, 0.1876589 , 0.23085264, 0.28398834,
       0.34935437, 0.42976579, 0.52868563, 0.650374  , 0.80007158,
       0.98422526, 1.21076588, 1.4894497 , 1.8322786 , 2.25401694])

and are given in units of :math:`0.05~r_{500}`, so that ``radial_bin_centers`` is dimensionless.

This array is used in the radial profiles plots (e.g. Fig. 4) as the :math:`x`-axis quantity. For
all profiles, gas particles are selected above :math:`10^5` K.

Density profiles
----------------

The density profiles are computed as the ratio between the hot gas mass in each shell and the
volume of the shell, and given in units of :math:`g/cm^3`. The ``.npy`` files contain a 2D array
with 30 spherical shells (for each profile) for all 200 snapshots.

>>> density_profiles = np.load("data/VR2915_+1res_ref/density_profiles.npy")
>>> density_profiles.shape
(200, 30)
>>> # The last row corresponds to the z=0 snapshot
>>> density_profiles[-1]


Temperature profiles
--------------------

Temperature profiles are mass-weighted, and given in Kelvin.

>>> temperature_profiles = np.load("data/VR2915_+1res_ref/temperature_profiles.npy")  # In Kelvin

To find the *scaled* (or self-similar) temperature profile, this array can be divided by the
virial (self-similar) temperature, :math:`T_{500}`:

>>> t500 = np.load("data/VR2915_+1res_ref/t500.npy")  # In Kelvin
>>> temperature_profiles_scaled = temperature_profiles / t500  # Dimensionless


Entropy profiles
----------------

Entropy profiles are computed from density and (mass-weighted) temperature profiles. For the
:math:`i^{th}` spherical shell, the entropy :math:`K_i` is given by

.. math::

    K_i = k_B T_i / n_{e,i}^{2/3},

where :math:`n_e` is the electron number density, which can be computed from the density profiles.
In this work, we assumed fully ionised gas to compute :math:`\rho \longrightarrow n_e`.

>>> entropy_profiles = np.load("data/VR2915_+1res_ref/entropy_profiles.npy")  # In keV cm^2

To find the *scaled* (or self-similar) entropy profile, this array can be divided by the
virial (self-similar) entropy, :math:`K_{500}`:

>>> k500 = np.load("data/VR2915_+1res_ref/k500.npy")  # In keV cm^2
>>> entropy_profiles_scaled = entropy_profiles / k500  # Dimensionless