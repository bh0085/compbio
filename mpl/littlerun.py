import mpl.littlegui as lg
import pickle

pressfun = pickle.load(open('pressfun.pickle'))
aw = lg.getWindow(pressfun)
aw.show()
