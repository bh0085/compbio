def ensemble2(data = None):
    import orange, orngTree, orngEnsemble

    if not data:
        data = orange.ExampleTable('bupa.tab')

    forest = orngEnsemble.RandomForestLearner(trees=50, name="forest")
    tree = orngTree.TreeLearner(minExamples=2, mForPrunning=2, \
                            sameMajorityPruning=True, name='tree')
    learners = [tree, forest]

    import orngTest, orngStat
    #results = orngTest.leaveOneOut(learners, data)
    results = orngTest.crossValidation(learners, data, folds=2)
    print "Learner  CA     Brier  AUC"
    for i in range(len(learners)):
        print "%-8s %5.3f  %5.3f  %5.3f" % (learners[i].name, \
                                                orngStat.CA(results)[i],
                                            
                                            orngStat.AUC(results)[i])
