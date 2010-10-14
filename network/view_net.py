#A sample QT GUI within which to view matrices.

import netutils as nu
import netsvd as ns
import matplotlib.pyplot as plt
from numpy import *

import sys, os, random
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class Fig(FigureCanvas):
    def __init__(self, parent = None, width = 5, height = 4, dpi = 100):
        fig = Figure(figsize = (width, height), dpi = dpi)
        self.axes = fig.add_subplot(111)
        self.axes.hold(False)
        
        #Init this canvas with the figure
        FigureCanvas.__init__(self,fig)
        self.setParent(parent)

class ViewApp(QMainWindow):
    def __init__(self,parent=None):
        QMainWindow.__init__(self,parent)
        
        self.data = arange(10.)

        self.create_figure()
        self.on_draw()
        print 'making app'

    def on_press(self,event):
        #box_points = event.artist.get_bbox().get_points()
        print dir(event)
        #print box_points

    def on_draw(self):
        """ Redraws the figure
        """
        self.fig.axes.clear()
        y = self.data
        x = range(len(y)) 
        self.fig.axes.plot(x,y)
        self.fig.draw

    def create_figure(self):
        #establish a frame
        self.main_frame = QWidget()
        #create the figure
        self.fig = Fig(self)
        #handle click events
        self.fig.mpl_connect('button_press_event', self.on_press)

        #overkill to add a layout at the moment but... why not?
        vbox = QVBoxLayout()
        vbox.addWidget(self.fig)        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
def view_gui(data):
    app = QApplication(sys.argv)
    print 'viewing!'
    form = ViewApp()
    form.data = data
    form.show()
    app.exec_()

def view():
    #arr = reshape(arange(100.),(10,10))
    arr = arange(100.)
    view_gui(arr)


def view_array(M):
    
    m = M.shape[0]
    n = M.shape[1]
    

    xs,ys,rs,cs = [], [], [], []
    max_r = 50^2
    rscl = max_r / (abs(M).max())
    
    
    for i in range(m):
        for j in range(n):
            y = i
            x = j
            e = M[i,j]
            r = abs(e) * rscl
            c = [int(e > 0), 0 , int(e < 0)]
            xs.append(x)
            ys.append(y)
            rs.append(r)
            cs.append(c)

    fig = plt.figure(1)
    fig.clear()
    sp = fig.add_subplot(221,frame_on = False); sp.set_axis_off()
    sp.scatter(xs,ys,rs,cs)
    plt.show()
