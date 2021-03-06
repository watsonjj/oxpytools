.. _getting_started:

*************************
Installation Instructions
*************************

This guide assumes you are using the *Anaconda* python distribution, installed locally. This handles the creation of virtual environments which will allow you to control the dependency libraries independently of other python projects. The creation and initiallisation of this virtual environment is covered in the installation instructions for ctapipe. A python 3.5 environment is required for this project.

------------
Dependencies
------------

* ctapipe: https://cta-observatory.github.io/ctapipe/getting_started/index.html
* targetio (python install): https://forge.in2p3.fr/projects/target/repository/targetio

Additional pip modules:

.. code-block:: bash

    pip install tqdm

--------------------------
Get the oxpytools software
--------------------------

Check out the `oxpytools <https://github.com/watsonjj/oxpytools>`__ repo:

.. code-block:: bash

    mkdir oxpytools
    git clone https://github.com/watsonjj/oxpytools.git

Now setup this checked out version for development:

.. code-block:: bash

    cd oxpytools
    make init     # will fetch required sub-repos and set up package
    make develop  # will make symlinks in your python library dir


Make sure the tests and examples code finds the test and example files.
Run the tests to make sure everything is OK:

.. code-block:: bash

   make test

Build the HTML docs locally and open them in your web browser:

.. code-block:: bash

   make doc-show

To update to the latest development version (merging in remote changes
to your local working copy):

.. code-block:: bash

   git pull

---------------
Developing Code
---------------

Checking out oxpytools in the manner described above is read-only, meaning that if you want to commit a change, you cannot (the master repo is locked to only the managers). Therefore, in order to develop, you need to make a personal fork on GitHub.
This is described in the AstroPy documentation http://astropy.readthedocs.org/en/latest/development/workflow/get_devel_version.html#get-devel .  You would need to of course change any reference to "astropy" the package to "oxpytools" and "astropy" the organization to "watsonjj", but the instructions should work.

Even easier (if you are on a Mac computer) is to use the `github-desktop GUI <https://desktop.github.com/>`_, which can do all of it for you automatically. It will handle the forking, syncing, and even allow you to issue pull-requests.
