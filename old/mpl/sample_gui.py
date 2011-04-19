def pressed(e):
    x = e.xdata
    y = e.ydata
    f = e.canvas.figure
    f.clear()
    a = f.add_axes([0,0,1,1],xlim = [0,10],ylim = [0,10])
    a.scatter(x,y,400,'black')
    e.canvas.draw()

def set_gui():
    from mpl import littlegui
    littlegui.setPressFun(pressed)
    
