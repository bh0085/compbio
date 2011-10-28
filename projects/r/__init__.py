"""
short demo.

"""

def plot_test():
   from rpy2.robjects.packages import importr
   graphics = importr('graphics')
   grdevices = importr('grDevices')
   base = importr('base')
   stats = importr('stats')
   
   import array
   
   x = array.array('i', range(10))
   y = stats.rnorm(10)
   
   grdevices.X11()
   
   graphics.par(mfrow = array.array('i', [2,2]))
   graphics.plot(x, y, ylab = "foo/bar", col = "red")
   
   kwargs = {'ylab':"foo/bar", 'type':"b", 'col':"blue", 'log':"x"}
   graphics.plot(x, y, **kwargs)
   
   
   m = base.matrix(stats.rnorm(100), ncol=5)
   pca = stats.princomp(m)
   graphics.plot(pca, main="Eigen values")
   stats.biplot(pca, main="biplot")

def mat_test():
   from rpy2.robjects import NA_Real
   from rpy2.rlike.container import TaggedList
   from rpy2.robjects.packages import importr
   
   base = importr('base')
   
   # create a numerical matrix of size 100x10 filled with NAs
   m = base.matrix(NA_Real, nrow=100, ncol=10)
   
   # fill the matrix
   for row_i in xrange(1, 100+1):
       for col_i in xrange(1, 10+1):
           m.rx[TaggedList((row_i, ), (col_i, ))] = row_i + col_i * 100


def bicluster():
    from rpy2.robjects import NA_Real
    from rpy2.rlike.container import TaggedList
    from rpy2.robjects.packages import importr
    
    base = importr('base')
    # create a numerical matrix of size 100x10 filled with NAs
    m = base.matrix(  #NA_Real, nrow=100, ncol=10)
    
    # fill the matrix
    for row_i in xrange(1, 100+1):
       for col_i in xrange(1, 10+1):
           m.rx[TaggedList((row_i, ), (col_i, ))] = row_i + col_i * 100

    
    biclust = importr("biclust")
    return biclust
