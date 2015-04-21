# Introduction #

Add your content here.

# System requirements #

In order to run **assemblyline** you must have the following installed:

  * A [python 2.7+](http://python.org) installation. We **strongly** recommend the [Enthought Python Distribution (EPD)](http://www.enthought.com/products/epd.php)
    * [numpy](http://numpy.scipy.org)
    * [networkx](http://networkx.lanl.gov)
  * [R](http://www.r-project.org)
    * `rpart` package (classification trees)
    * [ROCR](http://rocr.bioinf.mpi-sb.mpg.de) package (classifier visualization)

# Installation #

  * Download **assemblyline\_vXX.tar.gz** tarball from Downloads page
  * Uncompress tarball e.g. `tar -zxvf assemblyline_vXX.tar.gz`
```
$ tar -zxvf assemblyline_vXX.tar.gz
```
  * Change to the newly created directory
```
$ cd assemblyline_vXX
```
  * Run setup.py program to install assemblyline as follows
```
$ python setup.py build_ext --inplace
```
  * Add the assemblyline directory to your `PYTHONPATH` environment variable
```
$ export PYTHONPATH=<path_to_assemblyline>/assemblyline:$PYTHONPATH
```
  * Test that the **assemblyline** module can be imported in python
```
$ python
>>> import assemblyline
>>> 
```


# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages