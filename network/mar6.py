import compbio.network.netutils as nu
import numpy as np

trgs, tfs = nu.parse_net()
gname = 'FBgn0036754'
tg = trgs[gname]
weights = tg['weights']
wsort = np.argsort(weights)
tf_sel = [tg['tfs'][i] for i in wsort[:8]]

orders = [[1,5],
          [1,6],
          [7,0],
          0,2,3,4,5]

for o in orders[0:3]:
  ctfs = [tf_sel[i] for i in o]

  tg_intersection = [ tg 
                      for tg in tfs[ctfs[0]]['targets'] 
                      if tg in tfs[ctfs[1]]['targets']]
  tg_union = [tfs[ctfs[0]]['targets'] ]
  tg_union = tg_union + [ tg for tg in tfs[ctfs[1]]['targets'] if not tg in tg_union]

  print float(len(tg_intersection)) / float(len(tg_union))

o2  = [[[i,j] for i in range(8) if not i == j] for j in range(8)]
o2prime = []
for o in o2:
  o2prime.extend(o)
o2 = o2prime

similarities = []
for o in o2:
  ctfs = [tf_sel[i] for i in o]

  tg_intersection = [ tg 
                      for tg in tfs[ctfs[0]]['targets'] 
                      if tg in tfs[ctfs[1]]['targets']]
  tg_union = [tfs[ctfs[0]]['targets'] ]
  tg_union = tg_union + [ tg for tg in tfs[ctfs[1]]['targets'] if not tg in tg_union]
  similarities.append(float(len(tg_intersection)) / float(len(tg_union)))

mean_sim = np.mean(np.array(similarities))
print mean_sim


print len(tfs)
print tf_sel

