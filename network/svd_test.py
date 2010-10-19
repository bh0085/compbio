def svd_test():
    import mpl.littlegui as lg
    win = lg.makewindow(pressed)

def pressed():

    M = reshape(array([1,1,0,0,1,0,0,0,0][::-1]),(3,3))
    U, S, Vh = lin.svd(M)
    V = Vh.T
    
    

    f = plt.figure(0)
    f.clear()
    sp = f.add_subplot(221,frame_on = False);  sp.set_axis_off()
    sp.imshow(M , aspect = 'auto',interpolation = 'nearest')
    
    sp = f.add_subplot(222,frame_on = False);  sp.set_axis_off()
    sp.imshow(abs((U)) , aspect = 'auto',interpolation = 'nearest')
           
    sp = f.add_subplot(223,frame_on = False);  sp.set_axis_off()
    sp.imshow(abs(V) , aspect = 'auto',interpolation = 'nearest')
 
    sp = f.add_subplot(224,frame_on = False);  sp.set_axis_off()
    sp.imshow(diag(S), aspect = 'auto',interpolation = 'nearest')

    print V
    print U
    print S

    print dot(dot(U,diag(S)),V.T)

    print V
