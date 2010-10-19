import inspect
def parent():
    child()

def child():
    print inspect.stack()[0][3]
    print inspect.stack()[1][3]
