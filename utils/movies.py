import os, inspect
import re

def number_files(moviedir = os.getcwd()):
    for root, dirs, files in os.walk(moviedir):
        for f in files:
            if '.png' in f:
                f0 = os.path.join(root, f)
                num = re.search(re.compile( '\d+'), f)
                if num:
                    g = num.group()
                    newname = f.replace(g,'{0:08}'.format(int(g)))
                    os.rename(f0, os.path.join(root,newname))

def make_movie(moviedir = os.getcwd()):
  #
  # Now that we have graphed images of the dataset, we will stitch them
  # together using Mencoder to create a movie.  Each image will become
  # a single frame in the movie.
  #
  # We want to use Python to make what would normally be a command line
  # call to Mencoder.  Specifically, the command line call we want to
  # emulate is (without the initial '#'):
  # mencoder mf://*.png -mf type=png:w=400:h=300:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -o out.avi
  # See the MPlayer and Mencoder documentation for details.
  #

  import subprocess

  command = ('mencoder',
             'mf://{0}/*.png'.format(moviedir),
             '-mf',
             'type=png:w=400:h=400:fps=8',
             '-ovc',
             'lavc',
             '-lavcopts',
             'vcodec=mpeg4',
             '-oac',
             'copy',
             '-o',
             'out.avi')
  
  
  import subprocess
  #os.spawnvp(os.P_WAIT, 'mencoder', command)

  print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
  subprocess.check_call(command)

  print "\n\n The movie was written to 'output.avi'"
  
  print "\n\n You may want to delete *.png now.\n\n"
  

