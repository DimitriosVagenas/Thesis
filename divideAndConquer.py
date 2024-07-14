import algorithm as alg
import copy
import matplotlib.pyplot as plt
import sys

def constructPUV(xOfVertices,terrain,setP): #terrain = [[x1,x2,x3,...,xn],[y1,y2,y3,...,yn]]
    middle = int(len(terrain[0])/2)
    leftTerrain = [[],[]]
    rightTerrain = [[],[]]
    l = []
    r = []
    pL = [[],[]]
    plotInfo = []
    for i in range(len(terrain[0])):
        if middle >= i:
            leftTerrain[0].append(terrain[0][i])
            leftTerrain[1].append(terrain[1][i])
        if middle <= i:
            rightTerrain[0].append(terrain[0][i])
            rightTerrain[1].append(terrain[1][i])
    
    if len(leftTerrain[0])>6:
        l = constructPUV(xOfVertices, leftTerrain,setP)
        
    else:
        l = alg.pairPV(xOfVertices,copy.deepcopy(leftTerrain), setP, [])
    if len(rightTerrain[0])>6:
        r = constructPUV(xOfVertices, rightTerrain,setP)
    else:
        r = alg.pairPV(xOfVertices,copy.deepcopy(rightTerrain), setP, []) 
        
    
    vL = copy.deepcopy(leftTerrain)
    notPairPRVR = copy.deepcopy(rightTerrain)
    if notPairPRVR[0][0] == rightTerrain[0][0]:
        notPairPRVR[0].pop(0)
        notPairPRVR[1].pop(0)
    tempTerrain = copy.deepcopy(terrain) #current terrain without p points
    for pair in r:
        index = notPairPRVR[0].index(pair[1][0])
        notPairPRVR[0].pop(index)
        notPairPRVR[1].pop(index)
    for p in setP:
        if p[0] in notPairPRVR[0]:
            index = notPairPRVR[0].index(p[0])
            notPairPRVR[0].pop(index)
            notPairPRVR[1].pop(index)
        if p[0]>=leftTerrain[0][0] and p[0]<= leftTerrain[0][-1]:
            pL[0].append(p[0])
            pL[1].append(p[1])

    
    if len(pL[0])>0 and len(notPairPRVR[0])>0 :
        #a einai h koryfh v poy epitrepei sto pL shmeio na dei thn katakoryfh l
        a = alg.calculatePA(copy.deepcopy(pL),copy.deepcopy(vL)) ###########kenh lista
        
        #c einai to antistoixo a omws gia tis notPairPRVR koryfes
        c = alg.calculateVC(copy.deepcopy(rightTerrain),copy.deepcopy(notPairPRVR))
        interPA = [[],[]] #to shmeio pou temnei h eu8eia tou eu8ugrammou tmhmatos pL-a sthn l
        interVC = [[],[]] #to shmeio pou temnei h eu8eia tou eu8ugrammou tmhmatos VR'-c sthn l
        if len(a[0])>0 and len(c[0])>0:
            for i in range (len(a[0])):
                l1 = [[pL[0][i], pL[1][i]],[a[0][i],a[1][i]]]
                l2 = [[rightTerrain[0][0], 0], [rightTerrain[0][0] ,10]]
                temp = alg.intersectionPoint(l1,l2)
                
                interPA[0].append(temp[0])
                interPA[1].append(temp[1])
            for j in range (len(c[0])):
                l1 = [[notPairPRVR[0][j],notPairPRVR[1][j]],[c[0][j],c[1][j]]]
                l2 = [[rightTerrain[0][0], 0], [rightTerrain[0][0] ,10]]
                temp = alg.intersectionPoint(l1,l2)
                interVC[0].append(temp[0])
                interVC[1].append(temp[1])
            interPA[0] = [x-rightTerrain[0][0] for x in interPA[0]]
            interVC[0] = [x-rightTerrain[0][0] for x in interVC[0]]
            pL[0] = [x-rightTerrain[0][0] for x in pL[0]]
            notPairPRVR[0] = [x-rightTerrain[0][0] for x in notPairPRVR[0]]

            wp = [[],[],[],[]]
            zq = [[],[],[],[]]
            for i in range (len(interPA[0])):
                mnwp = alg.calculateLineFrom2Points([interPA[0][i],interPA[1][i]],[pL[0][i],pL[1][i]])
                wp[0].append(mnwp[0])
                wp[1].append(mnwp[1])
                #ray gp is y = -pL[0][i]*x + pL[1][i] for x >= m
                #starting point of wp = mnwp
            
            for i in range (len(interVC[0])):
                if interVC[0][i]!=notPairPRVR[0][i]:
                    mnzq = alg.calculateLineFrom2Points([interVC[0][i],interVC[1][i]],[notPairPRVR[0][i],notPairPRVR[1][i]])
                    zq[0].append(mnzq[0])
                    zq[1].append(mnzq[1])

                else:
                    continue
                    #here vertice is on l so slope is infinite
            #calculate the infinity x of wp and zq
        
            wpMinX = min(wp[0])-1.0012
            zqMaxX = max(zq[0])+1.0012
            for i in range (len(interPA[0])):
                endX = zqMaxX
                endY = -pL[0][i]*endX + pL[1][i]
                wp[2].append(endX)
                wp[3].append(endY)
            for i in range (len(interVC[0])):
                endX = wpMinX
                endY = -notPairPRVR[0][i]*endX + notPairPRVR[1][i]
                zq[2].append(endX)
                zq[3].append(endY)
                
            pairs = []

            for i in range (len(zq[0])):
                pairPLVR =[-1,-1] 
                mincrossy = sys.float_info.max
                
                for j in range (len(wp[0])):
                    l1 = [[zq[0][i],zq[1][i]],[zq[2][i],zq[3][i]]]
                    l2 = [[wp[0][j],wp[1][j]],[wp[2][j],wp[3][j]]]
        
                    cross = alg.intersectionPoint(l1,l2)
                    
                    if cross[0]>= wp[0][j] and cross[0]<= zq[0][i]:
                        if abs(cross[1]-zq[1][i])<mincrossy:
                            mincrossy = abs(cross[1]-zq[1][i])
                            pairPLVR = [j,i]
        
                if pairPLVR[0] != -1:
                    pairs.append([[pL[0][pairPLVR[0]]+rightTerrain[0][0],pL[1][pairPLVR[0]]],[notPairPRVR[0][pairPLVR[1]]+rightTerrain[0][0],notPairPRVR[1][pairPLVR[1]]]])
            interPA[0] = [x+rightTerrain[0][0] for x in interPA[0]]
            interVC[0] = [x+rightTerrain[0][0] for x in interVC[0]]
            pL[0] = [x+rightTerrain[0][0] for x in pL[0]]
            notPairPRVR[0] = [x+rightTerrain[0][0] for x in notPairPRVR[0]]
            plotInfo.extend([terrain,interPA,interVC,pL,notPairPRVR,wp,zq])
            #plotDNQ(plotInfo)
            return l+r+pairs

    return l+r




def plotDNQ(plotInfo):
    #terrain[0],interPA[1],interVC[2],pL[3],notPairPRVR[4],wp[5],zq[6]
    #plotInfo[i][0]:x,plotInfo[i][1]:y
    #wp[0]:startx wp[1]:starty wp[2]:endx wp[3]:endy
    #zq ''
    plt.figure()
    plt.subplot(211)
    plt.plot(plotInfo[0][0],plotInfo[0][1], color = "green")
    maxYmid = max([max(plotInfo[1][1]),max(plotInfo[2][1]),max(plotInfo[0][1])])
    minYmid = min([min(plotInfo[1][1]),min(plotInfo[2][1]),min(plotInfo[0][1])])
    plt.plot([plotInfo[2][0][0], plotInfo[2][0][0]],[maxYmid, minYmid],linestyle = 'dotted',color = 'pink')
    for i in range (len(plotInfo[1][0])):
        plt.plot([plotInfo[1][0][i]], plotInfo[1][1][i], marker = ".", color = "red")
        plt.plot(plotInfo[3][0][i], plotInfo[3][1][i], marker = ".", color = "red")
        plt.plot([plotInfo[3][0][i],plotInfo[1][0][i]],[plotInfo[3][1][i], plotInfo[1][1][i]],linestyle = 'dotted',color = '0')
    for j in range (len(plotInfo[2][0])):
        plt.plot([plotInfo[2][0][j]], [plotInfo[2][1][j]], marker = ".", color = "blue")
        plt.plot(plotInfo[4][0][j], plotInfo[4][1][j], marker = ".", color = "blue")
        plt.plot([plotInfo[4][0][j],plotInfo[2][0][j]],[plotInfo[4][1][j], plotInfo[2][1][j]],linestyle = 'dotted',color = '0')
    plt.subplot(212)
    for i in range (len(plotInfo[5])):
        plt.plot([plotInfo[5][0]], [plotInfo[5][1]], marker = ".", color = "red")
        plt.plot([plotInfo[5][0],plotInfo[5][2]],[plotInfo[5][1],plotInfo[5][3]],linestyle = 'dotted',color = 'red')
    for i in range (len(plotInfo[6])):
        plt.plot([plotInfo[6][0]], [plotInfo[6][1]], marker = ".", color = "blue")
        plt.plot([plotInfo[6][0],plotInfo[6][2]],[plotInfo[6][1],plotInfo[6][3]],linestyle = 'dotted',color = 'blue')
    plt.show()