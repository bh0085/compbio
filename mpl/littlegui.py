

from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MyWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        #self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        self.main_widget = QtGui.QWidget(self)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
def pressed(e):
    print dir(e)
    print e.canvas
    e.canvas.figure.clear()
    ax = e.canvas.figure.add_axes([0,0,1,1],xlim = [0,10],ylim=[0,10])
    ax.scatter(e.xdata,e.ydata,50,'black')
    e.canvas.draw()

def getWindow(pressfun = pressed):
    
    aw = MyWindow()
    vbox = QtGui.QVBoxLayout()
    figc = FigureCanvas(Figure())
    vbox.addWidget(figc)
    f = figc.figure
    ax = f.add_axes([0,0,1,1],xlim = [0,10],ylim=[0,10])
    ax.plot(range(10))
    aw.main_widget.setLayout(vbox)
    print 'setting pressfun'
    figc.mpl_connect('button_press_event',pressfun)
    return aw

def setPressFun(function):
    import pickle
    pickle.dump(function,open('pressfun.pickle','w'))

#joyride.pfinance.com
#analyzethe.us
