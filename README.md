# ZEMAS: ZoNexus Electron Microscopy Acquisition Suite

![](https://github.com/tobyfrancis/ZEMAS/blob/main/assets/ZoNexus_logo.jpeg) 

ZEMAS, the ZoNexus Electron Microscopy Accquisition Suite, is an application for automated control and acquisition of various transmission electron microscopy (TEM) modalities using ZoNexus stages.
 
## Installation

After cloning the repository, one can use the package manager [pip](https://pip.pypa.io/en/stable/) to install the dependencies by entering the folder ("cd ZEMAS") and running "pip install -e .", after which the application can be opened using either QT Creator (suggested) or, inside the folder, "python zemas.py". 

Warning: on certain systems, one of [pymatgen](https://github.com/materialsproject/pymatgen)'s dependencies, [spglib](https://github.com/spglib/spglib) will not successfully install. Currently, ZEMAS does not have any dependency on spglib, so this error can be ignored. (Note to developers: the included "crystals" library has certain dependencies on spglib but these portions of the code are not being  used).

## License
[GPL3](https://www.gnu.org/licenses/gpl-3.0.en.html)

