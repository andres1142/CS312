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
import copy
import numpy as np
from TSPClasses import *
import heapq
import itertools


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution,
        time spent to find solution, number of permutations tried during search, the
        solution found, and three null values for fields not used for this
        algorithm</returns>
    '''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
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

    def greedy(self, time_allowance=60.0):
        pass

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints:
        max queue size, total number of states created, and number of pruned states.</returns>
    '''

    def branchAndBound(self, time_allowance=1000000.0):
        #  matrix = self.createMatrix()
        matrix = np.array([
            [np.inf, 9, np.inf, 8, np.inf],
            [np.inf, np.inf, 4, np.inf, 2],
            [np.inf, 3, np.inf, 4, np.inf],
            [np.inf, 6, 7, np.inf, 12],
            [1, np.inf, np.inf, 10, np.inf]
        ])
        lower_matrix = self.calculateLowerCostMatrix(matrix)
        #  bssf = self.greedyBSSF(matrix)
        bssf = 31
        pq = []
        state = 1

        heapq.heappush(pq, State(lower_matrix["lowerMatrix"], lower_matrix["lowerBound"], "S{}".format(state)))
        state += 1

        start_time = time.time()
        while len(pq) != 0:
            # Get the state with the lowest lower bound
            current_state = heapq.heappop(pq)
            if current_state.lower_cost < bssf:
                # Expand the current state into smaller states
                subProblems = self.createSubproblems(current_state, state)
                state = subProblems["state"]
                for sub_problem in subProblems["subProblems"]:
                    if sub_problem.lower_cost < bssf:
                        sub_problem.citiesFrom.append(sub_problem.citiesTo[-1])
                        heapq.heappush(pq, sub_problem)
                if len(sub_problem.citiesFrom) == 5 and sub_problem.lower_cost < bssf:
                    bssf = sub_problem.lower_cost
        return bssf

    def createSubproblems(self, current_state, state_id):
        sub_problems = []
        last_city_visited = current_state.citiesFrom[-1]

        for i in range(len(current_state.matrix[last_city_visited])):
            if not i in current_state.citiesFrom:
                sub_state = copy.deepcopy(current_state)
                sub_state.state_id = "S{}".format(state_id)
                sub_state.citiesTo.append(i)  # Add the city to the list of cities
                sub_state.matrix[i][last_city_visited] = np.inf  # Set the city to infinite
                sub_state.lower_cost += sub_state.matrix[last_city_visited][
                    i]  # Add the cost to the lower cost matrix of sub_state

                # Convert the current row and the current column to infinity of the sub_state
                sub_state.matrix[last_city_visited, :] = np.inf
                sub_state.matrix[:, i] = np.inf

                # Calculate the minimum value for each row and add it to the lower bound
                for row in range(len(sub_state.matrix)):
                    if not row in sub_state.citiesFrom:
                        #  Find the minimum value in the current row
                        min_row_val = np.min(sub_state.matrix[row])
                        sub_state.lower_cost += min_row_val

                        for column in range(len(sub_state.matrix[row])):
                            if row != column:
                                sub_state.matrix[row][column] = sub_state.matrix[row][
                                                                    column] - min_row_val if min_row_val != np.inf else np.inf

                # Calculate the minimum value for each column and add it to the lower bound
                for column in range(len(sub_state.matrix)):
                    if not column in sub_state.citiesTo:
                        #  Find the minimum value in the current column
                        min_column_val = np.min(sub_state.matrix[:, column])
                        sub_state.lower_cost += min_column_val

                        for row in range(len(sub_state.matrix)):
                            if row != column:
                                sub_state.matrix[row][column] = sub_state.matrix[row][
                                                                    column] - min_column_val if min_column_val != np.inf else np.inf
                sub_problems.append(sub_state)
                state_id += 1
        return {"subProblems": sub_problems, "state": state_id}

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number of solutions found during search, the
        best solution found.  You may use the other three field however you like.
        algorithm</returns>
    '''

    def fancy(self, time_allowance=60.0):
        pass

    '''
    This function creates a matrix of the distances between all cities
    '''

    def createMatrix(self):
        matrix = np.full((len(self._scenario.getCities()), len(self._scenario.getCities())), np.inf)
        for i in range(len(self._scenario.getCities())):
            for j in range(len(self._scenario.getCities())):
                if i != j:
                    matrix[i][j] = self._scenario.getCities()[i].costTo(self._scenario.getCities()[j])

        return matrix

    ''' 
    This function calculates the lower bound for the branch and bound algorithm
    '''

    def calculateLowerCostMatrix(self, matrix):
        lower_matrix = np.full((matrix.shape[0], matrix.shape[1]), np.inf)
        lower_bound = 0

        # Calculate the minimum value for each row and add it to the lower bound
        for i in range(len(lower_matrix)):
            # Find the minimum value in the current row
            min_row_val = np.min(matrix[i])
            lower_bound += min_row_val

            for j in range(len(lower_matrix[i])):
                if i != j:
                    lower_matrix[i][j] = matrix[i][j] - min_row_val

        # Calculate the minimum value for each column and add it to the lower bound
        for i in range(len(lower_matrix[0])):
            # Find the minimum value in the current column
            min_col_val = np.min(lower_matrix[:, i])
            lower_bound += min_col_val

            for j in range(len(lower_matrix)):
                if i != j:
                    lower_matrix[j][i] = lower_matrix[j][i] - min_col_val

        return {"lowerMatrix": lower_matrix, "lowerBound": lower_bound}

    '''
    This function calculates the BSSF cost using a greedy algorithm. 
    It is used to calculate the upper bound for the branch and bound algorithm
    '''

    def greedyBSSF(self, matrix):
        cities = list(range(len(matrix)))
        current_city = cities.pop(0)
        bssf_cost = float('inf')
        cost = 0
        count = 0

        while cities:
            next_city = min(cities, key=lambda city: matrix[current_city][city])
            cost += matrix[current_city][next_city]
            cities.remove(next_city)
            current_city = next_city
            count += 1

            if not cities and cost < bssf_cost:
                bssf_cost = cost
        return bssf_cost
