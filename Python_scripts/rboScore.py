#    target = ['q']
#    predicted = ['a', 'q', 'd', 'c']

# getRBOScore(target, predicted)
# line =  'GSK3B_KNOCKDOWN_206_GDS4305_DN,MAPK1,CSNK2A1,RPS6KA4,RPS6KA1,MAPK14,CDK2,TAF1,MAPK8,MAPK3,GSK3B,MAP2K1,CSNK2A2,CDK1,MAPKAPK2,MAP3K8,MAPK9,PRKDC,HIPK2,RPS6KA5,RPS6KA2\n',

def getRBOScore(predicted, target):
    # RBO score
    # From:
    def rboScore(l1, l2, p=0.98):
        """
            Calculates Ranked Biased Overlap (RBO) score.
            l1 -- Ranked List 1
            l2 -- Ranked List 2
        """
        if l1 == None: l1 = []
        if l2 == None: l2 = []
        sl, ll = sorted([(len(l1), l1), (len(l2), l2)])
        s, S = sl
        l, L = ll
        if s == 0: return 0
        # Calculate the overlaps at ranks 1 through l
        # (the longer of the two lists)
        ss = set([])  # contains elements from the smaller list till depth i
        ls = set([])  # contains elements from the longer list till depth i
        x_d = {0: 0}
        sum1 = 0.0
        for i in range(l):
            x = L[i]
            y = S[i] if i < s else None
            d = i + 1
            # if two elements are same then
            # we don't need to add to either of the set
            if x == y:
                x_d[d] = x_d[d - 1] + 1.0
            # else add items to respective list
            # and calculate overlap
            else:
                ls.add(x)
                if y != None: ss.add(y)
                x_d[d] = x_d[d - 1] + (1.0 if x in ss else 0.0) + (1.0 if y in ls else 0.0)
                # calculate average overlap
            sum1 += x_d[d] / d * pow(p, d)
        sum2 = 0.0
        for i in range(l - s):
            d = s + i + 1
            sum2 += x_d[d] * (d - s) / (d * s) * pow(p, d)
        sum3 = ((x_d[l] - x_d[s]) / l + x_d[s] / s) * pow(p, l)
        # Equation 32
        rbo_ext = (1 - p) / p * (sum1 + sum2) + sum3
        return rbo_ext

    # Get a P that scales to the length of the predicted list
    import numpy as np
    def getP(d):

        w = 0.0
        p = 0.999
        print("Starting at W="+str(w)+"; P="+str(p))
        w_list=[]; p_list=[]
        print()
        while (w < 0.90 and p>0.01): # Adjust max w to change the % of weight that's captured in up until the depth d
            sumX = 0
            for i in np.arange(1, d):
                sumX += ((p ** i) / i)
            w = 1 - (p ** (d - 1)) + ((1 - p) / p) * d * (np.log(1 / (1 - p)) - sumX)
            p -= 0.01
            w_list.append(w)
            p_list.append(p)
            # if w!=w_list[-1]:
            #     break
        return (round(p_list[-2], 2))

    # Iterate over all possible combinations of target list and take the one that gives the best RBO score
    import itertools
    lst=[]
    P = getP(len(predicted))
    for itterTarget in list(itertools.permutations(target)):
        print("test 1")
        lst.append( rboScore(predicted, list(itterTarget), P) )
        print("test 2")
    print("P="+str(P)+";  RBOscore="+str(max(lst)))
    return max(lst)
