.. _properties:

Properties
=====

Data serialisation
------------

Object properties are stored in ``.npy`` files and serialised with NumPy. The following I/O
functions can be used:

.. code-block:: python

   import numpy as np

   my_array = np.ones(10)
   np.save(my_array, "my_array.npy")

   # Load it from file
   my_array = np.load(my_array)


Example: halo masses
----------------

For example:

>>> import numpy as np
>>> m200_z0 = np.load("<path-to>/data/VR18_+1res_ref/final_m200_0199.npy")
>>> print(m200_z0)
['shells', 'gorgonzola', 'parsley']

For example:

>>> import numpy as np
>>> m500 = np.load("<path-to>/data/VR18_+1res_ref/m500.npy")
>>> print(m500)