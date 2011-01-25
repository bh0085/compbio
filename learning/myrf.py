import  milk.supervised.randomforest as rf
import  milk.supervised.randomforest as rf
import milk

import orngTest
import orange
import orngTree, orngEnsemble
import orngTest, orngStat


import numpy as np
from numpy import *

class OLearn():
    def __init__(self,name, train_examples, measure = "retis", test_ex = None):
        self.training = train_examples
        self.name = name
        if name == 'knn':
            k = 5
            c = orange.kNNLearner(train_examples, k = k)
        elif name == 'tree':
            c = orange.TreeLearner(train_examples,measure = 'retis' ,
                                   mForPruning=4, minExamples=2)        
        
            if test_ex:
                print getpred(c, test_ex)
        elif name == 'forest':
            tree = orange.TreeLearner(measure=measure, 
                                      mForPruning=2, minExamples=4)

            c = orngEnsemble.RandomForestLearner(train_examples, 
                                                 trees=50, 
                                                 learner = tree)      
        
        self.classifier = c
    def predictions(self,ex):
        p = []
        
        print ex[1]
        #handle the case of ensemble learning
        if self.name in ['forest']:
            for t in ex:
                avg = array([c(t) for c in self.classfier.classifiers])
                avg = [float(a) for a in avg]
                p.append(np.mean(avg))
                p = array(p)
        else:
            p = getpred(self.classifier, ex)
        
        return p

        
def run_rf(x, y):
    

    x2 = array(x, float)
    y2 = array(y, int)

    #x2[x>0] = True
    #x2[x<=0] = False
    y2[y>0] = 1
    y2[y<=0] = 0

    features = squeeze(x2.T)
    labels = squeeze(y2.T)


    f2 = np.random.randn(100,20)
    f2[:50] *= 2
    l2 = np.repeat((0,1), 50)
    

    t = 'rf'
    if t == 'rf':
        obj = rf.rf_learner()
    else:
        obj = milk.supervised.tree.tree_learner()

    model = obj.train(features, labels)
    predictions = []
    for fs in features:
        predictions.append(  model.apply(fs))
    predictions = array(predictions)


    print array(labels,int)
    print array(predictions,int)


    #print shape(f2), shape(features)
    #print shape(l2), shape(labels)
    #print predictions
    #print shape(predictions)

    # raise Exception()

    return predictions

def printTree0(node, level):
    if not node:
        print " "*level + "<null node>"
        return

    if node.branchSelector:
        nodeDesc = node.branchSelector.classVar.name
        nodeCont = node.distribution
        print "\n" + "   "*level + "%s (%s)" % (nodeDesc, nodeCont),
        for i in range(len(node.branches)):
            print "\n" + "   "*level + ": %s" % node.branchDescriptions[i],
            printTree0(node.branches[i], level+1)
    else:
        nodeCont = node.distribution
        majorClass = node.nodeClassifier.defaultValue
        print "--> %s (%s) " % (majorClass, nodeCont),

def printTree(x):
    if isinstance(x, orange.TreeClassifier):
        printTree0(x.tree, 0)
    elif isinstance(x, orange.TreeNode):
        printTree0(x, 0)
    else:
        raise TypeError, "TreeClassifier or TreeNode expected"


import orange, orngTree


def get_ex(x,y,xb = False, yb = False):
    
    
    ctype = orange.FloatVariable()
    btype = orange.EnumVariable(values = orange.StringList(['0','1']))
    if xb: xtype = btype
    else: xtype = ctype
    if yb: ytype = btype
    else: ytype = ctype

    #note that ny is not set.
    #its gotta be 1
    nx = len(x)
    
    #map x, y onto binary characters if required
    if xb:
        xdat = array(x,str)
        xdat[less_equal(x,0)] = '0'
        xdat[greater(x,0)] = '1'
    else:
        xdat = array(x)
    if yb:
        ydat = array(y,str)
        ydat[less_equal(y,0)] = '0'
        ydat[greater(y,0)] = '1'
    else:
        ydat = array(y)
        
    data = [list(elt) for elt in xdat.T]
    
    if len(shape(y)) != 1:
        raise Exception("Sorry, y has to be 1 dimensional")
    nt = len(y)

    #stick ys on the end of the data example
    for i in range(nt):
        data[i].append(ydat[i])
    domain = orange.Domain([xtype]*nx, ytype )
    examples = orange.ExampleTable([orange.Example(domain, elt) for elt in data])

    return examples

def examples_from_inds(x,y,inds):
    all_examples = get_ex(x,y)
    train_examples = all_examples.getitems([int(x) for x in inds])    
    return train_examples

def getpred(c, ex):
    p = []
    for e in ex:
        p.append(c(e))
    return array(p,float)

def run_knn(x,y, 
                   train_idxs, 
                   test_idxs):
    all_examples = get_ex(x,y)
    train_examples = all_examples.getitems(train_idxs)
    test_examples = all_examples.getitems(test_idxs)
    
    k = 5
    c = orange.kNNLearner(train_examples, k = k)

    return getpred(c,test_examples)

def run_tree(x,y, 
             train_idxs, 
             test_idxs,
             measure = "retis"):
    all_examples = get_ex(x,y)
    train_examples = all_examples.getitems(train_idxs)
    test_examples = all_examples.getitems(test_idxs)
    
    raise Exception()

    k = 5
    c = orange.TreeLearner(train_examples,measure=measure, 
                           mForPruning=2, minExamples=4)

    return getpred(c,test_examples)

def run_forest(x,y, 
             train_idxs, 
             test_idxs,
             measure = "retis"):
    all_examples = get_ex(x,y, yb = True)
    train_examples = all_examples.getitems(train_idxs)
    test_examples = all_examples.getitems(test_idxs)
    
    k = 5
    tree = orange.TreeLearner(measure="retis", 
                              mForPruning=2, minExamples=4)

    forest = orngEnsemble.RandomForestLearner(train_examples, 
                                              trees=50, 
                                              learner = tree)
    p = []
    for t in test_examples:

        avg = array([c(t) for c in forest.classifiers])
        avg = [float(a) for a in avg]
        p.append(np.mean(avg))
    return p



def run(x,y):


    datab = get_ex(x,y,xb = 0, yb = 1)
    datac = get_ex(x,y,xb = 0, yb = 0)
    
    predictions = []


    classify = False
    regress = False
    ensemble = True
    if classify:
        data = datab
        lr = orange.LinearLearner(datab)
        tc = orange.TreeLearner(datab)
        lgr= orange.SVMLearner(datab)
        learners = [lr,tc,lgr]
    elif regress:
        data = datac
        maj = orange.MajorityLearner(train_data)
        maj.name = "default"
        rt = orngTree.TreeLearner(train_data, measure="retis", 
                                  mForPruning=2, minExamples=20)
        rt.name = "reg. tree"
        k = 5
        knn = orange.kNNLearner(train_data, k=k)
        knn.name = "k-NN (k=%i)" % k
        learners = [maj,rt, knn]

    elif ensemble:
        data = datab


        tree = orngTree.TreeLearner(mForPruning=2, name="tree")
        bs = orngEnsemble.BoostedLearner(tree, name="boosted tree")
        bg = orngEnsemble.BaggedLearner(tree, name="bagged tree")
        forest = orngEnsemble.RandomForestLearner(trees=50, name="forest")

        #data = orange.ExampleTable("lymphography.tab")

        learners = [tree, bs, bg,forest]
        results = orngTest.crossValidation(learners, data)
        print "Classification Accuracy:"
        for i in range(len(learners)):
            print ("%15s: %5.3f") % (learners[i].name, orngStat.CA(results)[i])

        
        raise Exception()
        tree = orngTree.TreeLearner( minExamples=2, mForPrunning=2, \
                                        sameMajorityPruning=True, name='tree')
        #qbs = orngEnsemble.BoostedLearner(tree,data, name="boosted tree")
        #bg = orngEnsemble.BaggedLearner(tree,data, name="bagged tree")
        forest = orngEnsemble.RandomForestLearner(data,trees=50, name="forest")





        #learners = [tree, bs, bg]
        #results = orngTest.crossValidation(learners, data)
        #print "Classification Accuracy:"
        #for i in range(len(learners)):
        #    print ("%15s: %5.3f") % (learners[i].name, orngStat.CA(results)[i])
        learners = [tree, forest]


        results = orngTest.crossValidation(learners, data)
        print "Learner  CA     Brier  AUC"
        for i in range(len(learners)):
            print "%-8s %5.3f  %5.3f  %5.3f" % (learners[i].name, \
                                                    orngStat.CA(results)[i],
                                               orngStat.BrierScore(results)[i],
                                                orngStat.AUC(results)[i])


    for l in learners:
        for d in data:
            p = l(d)
            
            if l == learners[-1]:
                predictions.append(p)

    predictions = array([predictions],float)
     
    print y
    print predictions
    return predictions
    
    

def run_otree(x, y):
    

    x2 = array(x, float)
    y2 = array(y, int)

    #x2[x>0] = True
    #x2[x<=0] = False
    y2[y>0] = 1
    y2[y<=0] = 0

    features = squeeze(x2.T)
    labels = squeeze(y2.T)

    data = orange.ExampleTable("lenses")
    treeClassifier = orange.TreeLearner(data)
    
    

    raise Exception()
    print array(labels,int)
    print array(predictions,int)


    return predictions
