#Documentation

To build the documentation pages, 
please have following packages installed within Python3

```
conda install sphinx
conda install sphinx_rtd_theme
conda install numpydoc
```

Then, in this folder `docs`, try:

```
make html
```

An HTTP format documentation is build in `build/html`, which you can open it by your web-browser. For example,

```
firefox build/html/index.html
```

