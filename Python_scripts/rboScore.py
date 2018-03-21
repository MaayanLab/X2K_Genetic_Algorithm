#    target = ['q']
#    predicted = ['a', 'q', 'd', 'c']

# target = targetKinases
# predicted = predictedKinases



# getRBOScore(target, predicted)
# line =  'GSK3B_KNOCKDOWN_206_GDS4305_DN,MAPK1,CSNK2A1,RPS6KA4,RPS6KA1,MAPK14,CDK2,TAF1,MAPK8,MAPK3,GSK3B,MAP2K1,CSNK2A2,CDK1,MAPKAPK2,MAP3K8,MAPK9,PRKDC,HIPK2,RPS6KA5,RPS6KA2\n'
# line = 'MAPK1_KNOCKDOWN_145_GSE31912_DN,MAPK14,MAPK1,MAPK3,CDK1,GSK3B,CDK2,MAPK9,RIPK2,PBK,PLK4,NLK,TSSK3,TRIM33,PRP4'
# target, predicted = parseTargetsPredicted(line, dataType)

def getRBOScore(predicted, target, staticD=False):
    # RBO score
    def rbo(l1, l2, p=0.98):
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
    def getP(d):
        import numpy as np
        w = 0.0
        p = 0.999
        #print("Starting at W=" + str(w) + "; P=" + str(p))
        #print("Evaluating P= "+str(p))
        while w < 0.90: # Adjust max w to change the % of weight that's captured in up until the depth d
            #print("Calculating new sumX...")
            sumX = 0
            for i in np.arange(1, d):
                sumX += ((p ** i) / i)
            #print("     sumX= "+str(sumX))
            #print("Calculating new w...")
            w = 1 - (p ** (d - 1)) + ((1 - p) / p) * d * (np.log(1 / (1 - p)) - sumX)
            #print("     w= "+str(w))
            p -= 0.001
        return (round(p, 2))

    # Iterate over all possible combinations of target list and take the one that gives the best RBO score
    import itertools
    lst=[]
    # Since all permutations will have the same length, can just calculate optimized P beforehand
    if staticD==False:
        P = getP(d=len(predicted)) # Dynamic
    else:
        P = getP(d=20) # Static
    # Repeat rbo for all permutations of the target list, and then simply take the highest scoring one
    ## ...because the target list is unranked, and therefore the order doesn't matter.
    for itterTarget in list(itertools.permutations(target)):
        #print("** itterTarget= " + str(itterTarget))
        lst.append( rbo(predicted, list(itterTarget), P) )
        #print("itterTarget: after append")
    rboScore = max(lst)
    #print("P="+str(P)+";  RBOscore="+str(rboScore))

    return rboScore
