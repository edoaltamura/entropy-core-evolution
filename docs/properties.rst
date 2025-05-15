.. _properties:

Properties
==========

The quantities describing global properties of the halos are provided in the ``data/`` directory.
Subdirectories are named as follows:

.. code-block:: text

    VR{halo-catalogue-id}_{resolution-level}_res_{subgrid-model}/

    # halo-catalogue-id: 2915 (group), or 18 (cluster)
    # resolution-level: -8 (low-res), or -1 (mid-res), or +8 (high-res, group only)
    # subgrid-model: ref (full-physics), or adiabatic (= non-radiative)

Let's take ``data/VR2915_+1res_ref/`` as an example for the code snippets in this guide.

Data serialisation
------------------

This simulation of the group with the Ref SWIFT-EAGLE model contains the following data arrays:

.. code-block:: bash

    $ ls ./

    baryon_mass.npy
    bh_edd_fraction.npy
    bh_edd_fraction_by_id.npy
    bh_mass.npy
    bh_mass_by_id.npy
    bh_mass_in_50kpc.npy
    centre_potential.npy
    ...


.. warning::

    Object properties are stored in ``.npy`` files and serialised with NumPy, which must be
    installed in your environment. The following I/O functions can be used:

    .. code-block:: python

       import numpy as np

       my_array = np.ones(10)
       np.save(my_array, "my_array.npy")

       # Load it from file
       my_array = np.load(my_array)


In the scripts below, we assume that:

- You have configured a Python environment

- The current working directory is the repository's root directory, ``entropy-core-evolution/``

- NumPy is imported:

>>> import numpy as np


Example: halo masses
--------------------

The self-similar masses are stored as arrays in files containing ``m200`` and ``m500`` in their
name. In general, the data is stored in a 1D array containing the VR catalogue halo mass for
every snapshot in the redshift output list (saved in ``redshifts.npy``). For convenience,
commonly used quantities like ``m200`` at :math:`z=0` are also duplicated in filenames starting
with ``final_`` and containing ``_0199``, which is the index of the final simulation snapshot
(out of 200).

>>> m200_z0 = np.load("data/VR2915_+1res_ref/final_m200_0199.npy")
>>> print(m200_z0)  # In solar masses
array(1.81625339e+13)


>>> m500 = np.load("data/VR2915_+1res_ref/m500.npy")
>>> print(m500.shape)  # In solar masses
(200,)
>>> print(m500[-1])
np.float64(8898277215097.906)

These data were used to produce Fig. 1 in the paper.

Example: other self-similar quantities
--------------------------------------

Similarly to :math:`M_{500}`, other self-similar quantities can be loaded from data files. A
list of available quantities can be found below.

- ``r500.npy``: the list of virial radii, in Mpc, for the main object at every snapshot.

- ``t500.npy``: the list of virial temperature, in K, for the main object at every snapshot.

- ``k500.npy``: the list of virial (self-similar) entropy, in :math:`kev~cm^2`, for the main
    object at every snapshot.


Example: baryon fractions
-------------------------

The hot (and cold) gas fractions are built from the halo mass, :math:`M_{500}`, and the hot (and
cold) gas mass computed by adding the mass of particles inside :math:`r_{500}` with a
temperature cut at :math:`10^5` K.

>>> m500 = np.load("data/VR2915_+1res_ref/m500.npy")  # In Solar masses
>>> m_hot = np.load("data/VR2915_+1res_ref/hot_gas_mass.npy")  # In Solar masses
>>> hot_gas_fraction = m_hot / m500
>>>
>>> # Return it every 20 snapshots
>>> print(hot_gas_fraction[::20])
array([0.00220512, 0.0002538 , 0.00587353, 0.00801783, 0.01054588,
       0.01692892, 0.02454117, 0.01514342, 0.03016408, 0.03807911])


Similarly, you can compute the other mass ratios used in Fig. 3 from the following files:

- ``cold_gas_mass.npy`` expresses the gas mass inside :math:`r_{500}` and below :math:`10^5` K.

- ``gas_mass.npy`` expresses the gas mass in Solar masses inside :math:`r_{500}` with any
    temperature.

- ``star_mass.npy`` expresses the total stellar mass, in Solar masses, inside :math:`r_{500}`.

- ``baryon_mass.npy`` expresses the total baryonic mass inside :math:`r_{500}`, obtained as the
    sum of the *total* gas and stellar masses above.

Finally, the central black hole mass in Solar masses is given in the ``bh_mass.npy`` file.
