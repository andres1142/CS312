#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


import time
import numpy as np
from TSPClasses import *
from heapQueue import heapQueue as queue
import itertools

def isInTimeLimit(start_time, time_allowance):
    return time.time()-start_time < time_allowance

class TSPSolver:
    def __init__( self, gui_view ):
        self._scenario = None

    def setupWithScenario( self, scenario ):
        self._scenario = scenario


    ''' <summary>
            This is the entry point for the default solver
            which just finds a valid random tour.  Note this could be used to find your
            initial BSSF.
        </summary>
        <returns>
            results dictionary for GUI that contains three ints: cost of solution,
            time spent to find solution, number of permutations tried during search, the
            solution found, and three null values for fields not used for this
            algorithm
        </returns>
    '''

    def defaultRandomTour( self, time_allowance=60.0 ):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and isInTimeLimit(start_time, time_allowance):
            # create a random permutation
            perm = np.random.permutation( ncities )
            route = []
            # Now build the route using the random permutation
            for i in range( ncities ):
                route.append( cities[ perm[i] ] )
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
            Takes the remaining unvisited cities and the route so far to find the next city greedily
            makes a priority queue with the cities with the distance from the last city as the key value
            takes the shortest path from that queue and adds it to the route while removing it from the list of unvisited cities
            total complexity is O(n*insert + deleteMin)
            with the heapQueue, that comes to O(n*log(n))
        </summary>
        <returns>
            [
                cities: the remaining unvisited cities
                route: the route so far
            ]
        </returns>
    '''
    def addNextGreedyCity(self, cities, route):
        q = queue(len(cities))
        # build queue (O(n*insert))
        for i in range(len(cities)):
            q.decreaseKey(i, route[-1].costTo(cities[i]))
        # add shortest option to the route (O(deleteMin))
        cityInd = q.deleteMin()
        route.append(cities[cityInd])
        del cities[cityInd]
        return [cities, route]

    ''' <summary>
        This is the entry point for the greedy solver, which you must implement for
        the group project (but it is probably a good idea to just do it for the branch-and
        bound project as a way to get your feet wet).  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number of solutions found, the best
        solution found, and three null values for fields not used for this
        algorithm</returns>
    '''

    def greedy( self,time_allowance=60.0 ):
        results = {}
        cities = self._scenario.getCities().copy()
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        route = [cities[0]]
        del cities[0]
        # complexity of O(n^2*log(n))
        while len(cities) > 0 and isInTimeLimit(start_time, time_allowance):
            [cities, route] = self.addNextGreedyCity(cities, route)
        bssf = TSPSolution(route)
        count += 1
        if bssf.cost < np.inf:
            # Found a valid route
            foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    def printMatrix(self, matrix):
        for arr in matrix:
            print(*arr, sep='\t')

    ''' <summary>
            reduce the given cost matrix
            complexity of O(n^2)
        </summary>
        <returns>
            [
                rcm: reduced cost matrix
                cost: the lower bound of the matrix
            ]
        </returns>
    '''
    def reduceMatrix(self, source):
        matrix = source.copy()
        cost = 0
        # reduce rows
        for i in range(len(matrix)):
            # find the minimum cost in the row
            minCost = np.inf
            unvisitedCount = 0
            for j in range(len(matrix[i])):
                if not np.isnan(matrix[i][j]):
                    unvisitedCount += 1
                    if matrix[i][j] < minCost:
                        minCost = matrix[i][j]
            if unvisitedCount > 0:
                if minCost == np.inf:
                    return [matrix, np.inf]
                cost = cost + minCost
            # remove the minimum from each row
            for j in range(len(matrix[i])):
                if not np.isnan(matrix[i][j]):
                    matrix[i][j] -= minCost

        # reduce columns
        for i in range(len(matrix[0])):
            # find the minimum cost in the column
            minCost = np.inf
            unvisitedCount = 0
            for j in range(len(matrix)):
                if not np.isnan(matrix[j][i]):
                    unvisitedCount += 1
                    if matrix[j][i] < minCost:
                        minCost = matrix[j][i]
            if unvisitedCount > 0:
                if minCost == np.inf:
                    return [matrix, np.inf]
                cost = cost + minCost
            # remove the minimum from each row
            for j in range(len(matrix)):
                if not np.isnan(matrix[j][i]):
                    matrix[j][i] -= minCost

        return [matrix, cost]

    ''' <summary>
            create the initial reduced cost matrix for an array of cities
            complexity of O(n^2)
        </summary>
        <returns>
            [
                rcm: reduced cost matrix
                cost: the lower bound of the matrix
            ]
        </returns>
    '''
    def initializeRCM(self, cities):
        rcm = np.array([[np.inf for _ in range(len(cities))] for _ in range(len(cities))])
        # O(n*n)
        for i in range(len(cities)):
            for j in range(len(cities)):
                rcm[i][j] = cities[i].costTo(cities[j])
        return self.reduceMatrix(rcm)


    ''' <summary>
            clears the row and column of the matrix
            complexity of O(n)
        </summary>
        <returns>
            the updated matrix
        </returns>
    '''
    def removeIJFromMatrix(self, source, i, j):
        matrix = source.copy()
        for k in range(len(matrix)):
            matrix[i][k] = None
            matrix[k][j] = None
        return matrix

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints:
        max queue size, total number of states created, and number of pruned states.</returns>
    '''

    def branchAndBound( self, time_allowance=60.0 ):
        results = {}
        cities = self._scenario.getCities().copy()
        ncities = len(cities)
        count = 0
        total = 0
        pruned = 0
        maxQueueSize = 0
        start_time = time.time()
        # start with the greedy solution as the BSSF
        # complexity of O(n^2*log(n))
        bssf = self.greedy()['soln']
        # find the initial RCM
        # complexity of O(n^2)
        [rcm, lowerBound] = self.initializeRCM(cities)
        matrices = [rcm]
        lowerBounds = [lowerBound]
        paths = [[0]]
        states = queue(0)
        states.insert(0,lowerBound)

        # this is potentially exponential, but we get a big speedup by the amount that we are able to prune 
        while not states.isEmpty() and isInTimeLimit(start_time, time_allowance):
            # choose a state O(deleteMin)
            # with heap, this is O(log(n))
            ind = states.deleteMin()
            # check that the state shouldn't be pruned
            if lowerBounds[ind] > bssf.cost:
                pruned += 1
                continue
            # expand the chosen state
            # this happens n times, so we get O(n^3) to expand
            node = paths[ind][-1]
            for i in range(ncities):
                # O(n^2) to copy
                matrix = matrices[ind].copy()
                if not np.isnan(matrix[node][i]) and matrix[node][i] < np.inf:
                    total += 1
                    lowerBound = lowerBounds[ind] + matrix[node][i]
                    newMatrix = self.removeIJFromMatrix(matrix, node, i)
                    # O(n^2) to reduce
                    [rcm, cost] = self.reduceMatrix(newMatrix)
                    lowerBound += cost
                    if lowerBound < bssf.cost:
                        key = len(matrices)
                        matrices.append(rcm)
                        lowerBounds.append(lowerBound)
                        path = paths[ind].copy()
                        path.append(i)
                        paths.append(path)
                        weight = len(path) * max( 10000 / (count + 1), 700 )
                        # add new state, balancing the lowerBound with the depth in the tree
                        # O(insert) depends on the queue
                        # because we're using a heap implementation, it is O(log(n))
                        states.insert(key, lowerBound - (weight * len(path)))

                        # check if we have a path
                        # O(n)
                        if len(paths[key]) == ncities:
                            count += 1
                            citiesOnPath = []
                            for ind in paths[key]:
                                citiesOnPath.append(cities[ind])

                            newSolution = TSPSolution(citiesOnPath)
                            if newSolution.cost < bssf.cost:
                                bssf = newSolution
                            else:
                                pruned += 1
                    elif lowerBound < np.inf:
                        pruned += 1
                maxQueueSize = max(maxQueueSize, states.size())

        end_time = time.time()
        # if time runs out, add all states with a lower bound greater than the BSSF
        if not isInTimeLimit(start_time, time_allowance):
            ind = states.deleteMin()
            while not states.isEmpty():
                if lowerBounds[ind] < bssf.cost: 
                    pruned += 1
                ind = states.deleteMin()
        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = maxQueueSize
        results['total'] = total
        results['pruned'] = pruned
        print(results)
        print(len(paths))
        return results

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number of solutions found during search, the
        best solution found.  You may use the other three field however you like.
        algorithm</returns>
    '''

    def fancy( self,time_allowance=60.0 ):
        pass
