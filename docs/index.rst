.. include:: ../README.rst

Data and scripts from the entropy-core evolution project
===================================

This documentation complements the GitHub repository with data and scripts by providing minimal
examples to help you **include data from Altamura et al. (2025) in your figures in no time!**

Check out the :doc:`science` section for further information on the results and implications for
the entropy-core problem. Visit the :ref:`properties`, :ref:`profiles`, :ref:`lagrangian-histories`, and
:ref:`projected-maps` sections of the documentation for minimal
Python examples to read in and plot the data from our paper directly into your figures. Each
section describes the serialisation of the data, the file structure, and provides abstract
I/O functions to read data from file.

Citing
--------

.. note::

    If using these results, data or scripts in your work, please consider citing Altamura et al.
    (2005) following the guidelines below.

Bibtex:

.. code-block:: bibtex

   bibtex

.. note::

    If you are also using results about the entropy-core problem or results at :math:`z=0` (equiv.
    local universe), please also cite
    `Altamura et al. (2023) <https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.3164A>`_ as follows.

Bibtex:

.. code-block:: bibtex

   @ARTICLE{2023MNRAS.520.3164A,
           author = {{Altamura}, Edoardo and {Kay}, Scott T. and {Bower}, Richard G. and {Schaller}, Matthieu and {Bah{\'e}}, Yannick M. and {Schaye}, Joop and {Borrow}, Josh and {Towler}, Imogen},
            title = "{EAGLE-like simulation models do not solve the entropy core problem in groups and clusters of galaxies}",
          journal = {\mnras},
         keywords = {hydrodynamics, methods: numerical, software: simulations, galaxies: clusters, galaxies: fundamental parameters, galaxies: groups - tions, Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies},
             year = 2023,
            month = apr,
           volume = {520},
           number = {2},
            pages = {3164-3186},
              doi = {10.1093/mnras/stad342},
    archivePrefix = {arXiv},
           eprint = {2210.09978},
     primaryClass = {astro-ph.CO},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.3164A},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

.. note::

    If using the *topological closure* technique to prepare contamination-free initial initial
    conditions for zoom-in cosmological simulations, please cite the following
    `Ph.D. thesis <https://ui.adsabs.harvard.edu/abs/2023PhDT.........8A>`_.

Bibtex:

.. code-block:: bibtex

   @PHDTHESIS{2023PhDT.........8A,
           author = {{Altamura}, Edoardo},
            title = "{Building models of the Universe with hydrodynamic simulations}",
         keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Astrophysics of Galaxies, Physics - Computational Physics, Physics - Fluid Dynamics},
           school = {University of Manchester, UK},
             year = 2023,
            month = dec,
           adsurl = {https://ui.adsabs.harvard.edu/abs/2023PhDT.........8A},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

Contents
--------

.. toctree::

   Home <self>
   science
   properties
   profiles
   lagrangian-histories
   projected-maps
