import inference

def getNet(name = 'easy0'):
    if name == 'easy0':
        return inference.get_easy0()
    else: raise Exception()
    return 

