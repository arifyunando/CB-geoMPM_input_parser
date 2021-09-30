from queue import Queue
import os


"""
Bresenham Algothms for drawing lines

use -> plotLine(pStart, pEnd)
p_start : <tuple> initial point
p_end   : <tuple> end point

return <list.tuple> in between points 
"""

def plotLineM(x0, y0, x1, y1):
    dx, dy = int(x1 - x0), int(y1 - y0)
    yi, dy = (-1, -dy) if dy < 0 else (1, dy)
    D  = 2*dy - dx
    y  = int(y0)

    rv = list()
    for x in range(int(x0), int(x1) + 1):
        rv.append((x, y))
        if D > 0:
            y += yi
            D += 2*(dy - dx)
        else:
            D += 2*dy
    return rv

def plotLineS(x0, y0, x1, y1):
    dx, dy = int(x1 - x0), int(y1 - y0)
    xi, dx = (-1, -dx) if dx < 0 else (1, dx)
    D  = 2*dx - dy
    x  = int(x0)

    rv = list()
    for y in range(int(y0), int(y1) + 1):
        rv.append((x, y))
        if D > 0:
            x += xi
            D += 2*(dx - dy)
        else:
            D += 2*dx
    return rv

def plotLine(pStart, pEnd):
    # Bresenham Algorithms for interger only domain
    x0, y0 = pStart
    x1, y1 = pEnd
    if abs(y1 - y0) < abs(x1 - x0):
        return plotLineM(x1, y1, x0, y0) if x0 > x1 else plotLineM(x0, y0, x1, y1)
    else:
        return plotLineS(x1, y1, x0, y0) if y0 > y1 else plotLineS(x0, y0, x1, y1)

"""
Fill Algorithm
use -> floodFill(sPoint, boundaryList)
"""

def domain(tuplist):
    x   = [i[0] for i in tuplist]
    y   = [i[1] for i in tuplist]
    rv  = []
    for i in range(min(y), max(y) + 1):
        for j in range(min(x), max(x) + 1):
           rv.append((j, i))
    domainDict = {k : 0 for k in rv}
    for i in tuplist:
        domainDict[i] += 1
    return domainDict

def fourWay(tup):
    return [(tup[0] + 1, tup[1]), (tup[0], tup[1] + 1), (tup[0] - 1, tup[1]), (tup[0], tup[1] -1)]
    
def floodFill(sPoint, domain):
    q_list = Queue()
    q_list.put(sPoint)
    
    while not q_list.empty():
        q = q_list.get()
        if domain[q] == 0:
            domain[q] += 1
            for direction in fourWay(q):
                q_list.put(direction)
        else:
            continue
    return [k for k, v in domain.items() if v > 0]

"""
mesh Class

+ set dimension
+ discretization scaling
+ output files
"""

class mesh():

    # Constructor
    def __init__(self, dx, dy = None, dz = None):
        self.delta = dx
        self.dx = dx
        self.dy = dx if dy == None else dy
        self.dz = dx if dz == None else dz


    # Private Methods
    def __fillNodes(self):
        domainDict  = domain(self.__defineBorder(self.__defineCorner()))
        centroid    = (int((self.x/self.dx)/2), int((self.y/self.dy)/2))
        nodes       = floodFill(centroid, domainDict)
        nodes.sort(key = lambda x : (x[1], x[0]))
        self.ints   = nodes[:]
        self.nodes  = nodes[:]
        for key, ele in enumerate(self.nodes):
            self.nodes[key] = [item*self.delta for item in ele]
        return self.nodes

    def __defineBorder(self, corners):
        border  = set()
        c       = corners[:] + corners[:1]
        for i in range(len(c) - 1):
            border.update(plotLine(c[i], c[i+1]))
        self.border = border
        return self.border

    def __defineCorner(self):
        rv = [
            (0, 0),
            (int(self.x/self.dx), 0),
            (int(self.x/self.dx), int(self.y/self.dy)),
            (0, int(self.y/self.dy))]
        return rv

    def __setCell(self):
        self.cell = []
        for node in self.ints:
            try:
                self.cell.append(
                    (
                        self.ints.index(node),
                        self.ints.index((node[0] + 1, node[1])),
                        self.ints.index((node[0] + 1, node[1] + 1)),
                        self.ints.index((node[0], node[1] + 1))
                    )
                )
            except:
                continue
        return self.cell


    # Public Methods
    def setMesh(self, x, y, z = 0):
        self.x, self.y, self.z = x, y, z
        self.__fillNodes()
        self.__setCell()

    def rotateMesh(self, tetha):
        pass

    def printFile(self, fileName = 'mesh', folderName = "project", digits = 0):
        path = f"result/{folderName}"
        try:
            os.mkdir(path)
        except OSError:
            print ("Path %s existed already, continue with writting files" % path)
        else:
            print ("Successfully created the directory %s " % path)
            
        with open(f'{path}/{fileName}.txt', 'w') as f:
            f.write(f'! {len(self.nodes)} nodes, {len(self.cell)} cells\n')
            f.write(f'\n{len(self.nodes)} \t{len(self.cell)}\n')
            # f.write('\n! Nodes\n')
            for i in self.nodes:
                f.write(' '.join([f'{j:.{digits}f}' for j in i]) + '\n')
                # f.write(''.join([f'{f"{j:.{digits}f}":>7}' for j in i]) + '\n')
            # f.write('\n! Cells\n')
            for i in self.cell:
                f.write(' '.join([f'{j}' for j in i]) + '\n')
                # f.write(''.join([f'{str(j):>7}' for j in i]) + '\n')

"""
particle Class

+ set cornerPoint
+ rotate
+ output files
"""

class particle():

    # Constructor
    def __init__(self, pointList, delta = 1, intPoint = None):
        self.delta = delta
        points = [(int(i/delta), int(j/delta)) for i, j  in pointList]
        self.__defineBorders(points)
        self.__fillNodes(intPoint)

    # Private Methods
    def __defineBorders(self, pointList):
        border  = set()
        p       = pointList[:] + pointList[:1]
        for i in range(len(p) - 1):
            border.update(plotLine(p[i], p[i+1]))
        self.border = border
        return self.border

    def __fillNodes(self, intPoint = None):
        domainDict  = domain(self.border)
        # print(domainDict)
        x           = [i[0] for i in self.border]
        y           = [i[1] for i in self.border]
        centroid    = (int(sum(x)/len(x)), int(sum(y)/len(y))) if intPoint == None else intPoint
        self.nodes  = floodFill(centroid, domainDict)
        self.nodes.sort(key = lambda x : (x[1], x[0]))
        for key, ele in enumerate(self.nodes):
            self.nodes[key] = [item*self.delta for item in ele]
        return self.nodes
    

    # Public Methods
    def setPoints(self, pointList):
        self.__defineBorders(pointList)

    def translate(self, x, y):
        displacements = (x, y)
        for key, ele in enumerate(self.nodes):
            self.nodes[key] = [item + d for item, d in zip(ele, displacements)]

    def rotate(self, tetha):
        pass

    def printFile(self, fileName = 'Particles', folderName = "project", digits = 0, spacing = 7):
        path = f"result/{folderName}"
        try:
            os.mkdir(path)
        except OSError:
            print ("Path %s existed already, continue with writting files" % path)
        else:
            print ("Successfully created the directory %s " % path)
            
        with open(f'{path}/{fileName}.txt', 'w') as f:
            f.write(f'! {len(self.nodes)} particle(s) \n')
            f.write(f'{len(self.nodes)} \n')
            for i in self.nodes:
                f.write(' '.join([f'{f"{j:.{digits}f}":>{spacing}}' for j in i]) + '\n')
    
    