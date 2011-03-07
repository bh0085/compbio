from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, \
     RotatingMarker, ReverseBar, SimpleProgress

def simple(maxval):
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=maxval).start()
    return pbar
