Installing ``whampy``
=====================

Requirements
------------

This package has the following dependencies:

* `Python <http://www.python.org>`_ 3.10 or later
* `Numpy <http://www.numpy.org>`_ 1.26 or later
* `Astropy <http://www.astropy.org>`_ 5.3 or later
* `SciPy <http://scipy.org>`_ 1.13 or later
* `PyNeb <https://github.com/Morisset/PyNeb_devel>`_ 1.1 or later
* `matplotlib <http://matplotlib.org>`_ 3.8 or later
* `seaborn <https://seaborn.pydata.org/index.html>`_ 0.12 or later

Installation
------------

You can install pyonized using pip::

  pip install pyonized

To stall the latest developer version of pyonized you can type::

    git clone https://github.com/Deech08/pyonized.git
    cd pyonized
    python setup.py install

You may need to add the ``--user`` option to the last line `if you do not
have root access <https://docs.python.org/2/install/#alternate-installation-the-user-scheme>`_.
You can also install the latest developer version in a single line with pip::

    pip install git+https://github.com/Deech08/pyonized.git
