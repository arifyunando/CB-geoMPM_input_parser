from queue import Queue
import threading
import os
import time

"""
Bresenham Algothms for drawing lines

use -> plotLine(pStart, pEnd)
p_start : <tuple> initial point
p_end   : <tuple> end point

return <list.tuple> in between points 
"""

def plotLineM(x0, y0, x1, y1):
    dx = int(x1 - x0)
    dy = int(y1 - y0)
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
    dx = int(x1 - x0)
    dy = int(y1 - y0)
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
    x = [i[0] for i in tuplist]
    y = [i[1] for i in tuplist]
    rv = []
    for i in range(min(y), max(y) + 1):
        for j in range(min(x), max(x) + 1):
           rv.append((j, i))
    domainDict = {k : 0 for k in rv}
    for i in tuplist:
        domainDict[i] += 1
    return domainDict


def fourWay(tup):
    return [(tup[0] + 1, tup[1]), (tup[0], tup[1] + 1), (tup[0] - 1, tup[1]), (tup[0], tup[1] -1)]
    
# def query(queue, domain):
#     if not queue.empty():    
#         q = queue.get()
#         if domain[q] == 0:
#             domain[q] += 1
#             for direction in fourWay(q):
#                 queue.put(direction)
#         else:
#             pass

def floodFill(sPoint, domain):
    q_list = Queue()
    q_list.put(sPoint)
    # cores = os.cpu_count()
    # condition = True
    while not q_list.empty():
        # processes = [threading.Thread(target=query, args=[q_list, domain]) for _ in range(cores)]
        # for process in processes:
        #     process.start()
        # for process in processes:
        #     process.join()
        q = q_list.get()
        if domain[q] == 0:
            domain[q] += 1
            for direction in fourWay(q):
                q_list.put(direction)
        else:
            continue
        # condition = not q_list.empty()

    return [k for k, v in domain.items() if v > 0]

"""
Sorting Tuples by its index
default index = 0
"""

def sortTuple(tup, index = 0):
    tup.sort(key = lambda x : x[index])
    return tup

"""
mesh Class

+ set dimension
+ discretization scaling
+ output files

"""

class mesh():
    """
    To use:
    1. Initialize
    2. setMesh -> define discretization size
    3. 
    """
    # Constructor
    def __init__(self, dx, dy = None, dz = None):
        self.dx = dx
        # self.dy = dx if dy == None else dy
        self.dz = dx if dz == None else dz

    # Private Methods
    def fillNodes(self):
        # start = time.perf_counter()
        domainDict = domain(self.defineBorder(self.defineCorner()))
        centroid = (int((self.x/self.dx)/2), int((self.z/self.dz)/2))
        # check1 = time.perf_counter()
        self.nodes = floodFill(centroid, domainDict)
        # print(f'check 1 : {check1 - start}\ncheck 2 : {time.perf_counter() - check1}')
        return self.nodes

    def defineBorder(self, corners):
        border = set()
        c = corners[:] + corners[:1]
        for i in range(len(c) - 1):
            border.update(plotLine(c[i], c[i+1]))
        self.border = border
        return self.border

    def defineCorner(self):
        rv = [0, 0, 0, 0]
        rv[0] = (0, 0)
        rv[1] = (int(self.x/self.dx), 0)
        rv[2] = (int(self.x/self.dx), int(self.z/self.dz))
        rv[3] = (0, int(self.z/self.dz))
        return rv

    def sortNodes(self):
        nodes = list(self.nodes)
        return sortTuple(sortTuple(nodes, 1), 0)


    # Public Methods
    def setMesh(self, x, z, y = 0):
        self.x, self.z = x, z
        # self.y = y

    def rotateMesh(self, tetha):
        pass

    def printFile(self, fileName = 'mesh.txt'):
        pass