<p align="center">
  <a href="https://github.com/jlsteenwyk/clipkit">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/ClipKIT/master/docs/_static/img/logo.jpg" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/ClipKIT/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/clipkit/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/clipkit/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/clipkit/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/clipkit">
        </a>
        <a href="https://github.com/JLSteenwyk/ClipKIT/actions" alt="Build">
            <img src="https://img.shields.io/github/workflow/status/jlsteenwyk/clipkit/CI/master">
        </a>
        <a href="https://codecov.io/gh/jlsteenwyk/clipkit" alt="Coverage">
          <img src="https://codecov.io/gh/jlsteenwyk/clipkit/branch/master/graph/badge.svg?token=0J49I6441V">
        </a>
        <a href="https://github.com/jlsteenwyk/clipkit/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/clipkit">
        </a>
        <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
            <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
                alt="follow on Twitter">
        </a>
    </p>
</p>

ClipKIT is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.<br /><br />
If you found clipkit useful, please cite *ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference*. Steenwyk et al. 2020, PLoS Biology. doi: [10.1371/journal.pbio.3001007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007).
<br /><br />

---

<br />

This documentation covers downloading and installing ClipKIT. Details about each function as well as tutorials for using ClipKIT are available in the [online documentation](https://jlsteenwyk.com/ClipKIT/).

<br />

**Installation**

**If you are having trouble installing ClipKIT, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [twitter](https://twitter.com/jlsteenwyk) to get help.**

To install using *pip*, we strongly recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:
```shell
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install clipkit
pip install clipkit
```
**Note, the virtual environment must be activated to use *clipkit*.**

After using ClipKIT, you may wish to deactivate your virtual environment and can do so using the following command:
```shell
# deactivate virtual environment
deactivate
```

<br />

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:
```shell
# download
git clone https://github.com/JLSteenwyk/ClipKIT.git
cd PhyKIT/
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install
make install
```
To deactivate your virtual environment, use the following command:
```shell
# deactivate virtual environment
deactivate
```
**Note, the virtual environment must be activated to use *clipkit*.**

<br />

To install via anaconda, execute the follwoing command:

``` shell
conda install -c jlsteenwyk clipkit
```
Visit here for more information: https://anaconda.org/jlsteenwyk/clipkit
