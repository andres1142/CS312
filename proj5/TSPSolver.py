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
import threading
import random


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
        results = {}
        cities = self._scenario.getCities()
        matrix = self.createMatrix()
        ncities = list(range(len(matrix)))
        current_city = ncities.pop(0)
        bssf_cost = float('inf')
        best_route = []
        cost = 0
        count = 0
        route = [cities[current_city]]

        start_time = time.time()
        while ncities and time.time() - start_time < time_allowance:
            next_city = min(ncities, key=lambda city: matrix[current_city][city])
            cost += matrix[current_city][next_city]
            ncities.remove(next_city)
            current_city = next_city
            count += 1
            route.append(cities[current_city])

            if not ncities and cost < bssf_cost:
                bssf_cost = cost
                best_route = route[:]
        end_time = time.time()

        results['cost'] = bssf_cost if bssf_cost != math.inf else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = TSPSolution(best_route)
        results['max'] = None
        results['total'] = None
        results['pruned'] = None

        return results

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints:
        max queue size, total number of states created, and number of pruned states.</returns>
    '''

    # TODO: Check how many childs I am creating thorughout the whole algorithm
    def branchAndBound(self, time_allowance=60.0):
        matrix = self.createMatrix()
        lower_matrix = self.calculateLowerCostMatrix(matrix)
        greedyResult = self.greedyBSSF(matrix)
        bssf = greedyResult["cost"]
        bssf_solution = greedyResult["solution"]
        pq = []
        state = 1
        max_queue_size = 0
        solutions = 0
        pruned = 0

        heapq.heappush(pq, State(lower_matrix["lowerMatrix"], lower_matrix["lowerBound"], "S{}".format(state)))
        state += 1

        start_time = time.time()
        while pq and time.time() - start_time < time_allowance:
            # Get the state with the lowest lower bound
            current_state = heapq.heappop(pq)
            if current_state.lower_cost < bssf:
                # Expand the current state into smaller states
                subProblems = self.createSubproblems(current_state, state)
                state = subProblems["state"]

                if len(current_state.citiesFrom) == len(self._scenario.getCities()) and current_state.lower_cost < bssf:
                    bssf_solution = TSPSolution(self.map_cities(current_state.citiesFrom))
                    bssf = current_state.lower_cost
                    solutions += 1
                for sub_problem in subProblems["subProblems"]:
                    if sub_problem.lower_cost < bssf:
                        heapq.heappush(pq, sub_problem)
                    else:
                        pruned += 1
                if len(pq) > max_queue_size:
                    max_queue_size = len(pq)
            else:
                pruned += 1

        end_time = time.time()
        result = {'cost': bssf, 'time': end_time - start_time, 'count': solutions, 'soln': bssf_solution,
                  'max': max_queue_size, 'total': state, 'pruned': pruned + len(pq)}
        return result

    def map_cities(self, citiesIndices):
        cities = self._scenario.getCities()
        city_map = []

        for i in range(len(citiesIndices)):
            city_map.append(cities[citiesIndices[i]])
        return city_map

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
                        min_row_val = np.min(sub_state.matrix[row])
                        if np.isfinite(min_row_val):
                            sub_state.lower_cost += min_row_val
                            sub_state.matrix[row, :] -= min_row_val  # subtract the minimum value from the entire row
                        else:
                            sub_state.lower_cost = np.inf
                            break

                # Calculate the minimum value for each column and add it to the lower bound
                for column in range(sub_state.matrix.shape[1]):
                    if not column in sub_state.citiesTo:
                        min_column_val = np.min(sub_state.matrix[:, column])
                        if np.isfinite(min_column_val):
                            sub_state.lower_cost += min_column_val
                            sub_state.matrix[:,
                            column] -= min_column_val  # subtract the minimum value from the entire column
                        else:
                            sub_state.lower_cost = np.inf
                            break

                sub_state.citiesFrom.append(sub_state.citiesTo[-1])
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
        return self.AntColonyOptimization(time_allowance)

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
        cities = self._scenario.getCities()
        ncities = list(range(len(matrix)))
        current_city = ncities.pop(0)
        bssf_cost = float('inf')
        cost = 0
        count = 0
        best_route = []
        route = [cities[current_city]]

        while ncities:
            next_city = min(ncities, key=lambda city: matrix[current_city][city])
            cost += matrix[current_city][next_city]
            ncities.remove(next_city)
            current_city = next_city
            count += 1
            route.append(cities[current_city])

            if not ncities and cost < bssf_cost:
                bssf_cost = cost
                best_route = route[:]
        return {"solution": TSPSolution(best_route), "cost": bssf_cost}

    def AntColonyOptimization(self, time_allowance=60.0):
        # Initialize counts
        results = {}
        count = 0
        start_time = time.time()

        # Initialize algorithm parameters
        num_ants = 50  # Number of ants per generation
        num_iterations = 300  # Number of iterations
        evaporation_rate = 0.095  # Evaporation rate
        alpha = 0.9  # Alpha: Pheromone influence
        beta = 1.5  # Beta: Heuristic influence
        generations_since_update = 0
        generations = 0

        # Create_distance_matrix
        matrix = self.createMatrix()
        # Initialize pheromone matrix with initial_pheromone values
        pheromone_matrix = np.ones((len(self._scenario.getCities()),
                                    len(self._scenario.getCities())))

        best_solution = None
        best_cost = np.inf

        while generations_since_update < num_iterations and time.time() - start_time < time_allowance:
            generations_since_update += 1
            generations += 1

            solutions = self.findAntsPaths(num_ants, pheromone_matrix, matrix, alpha, beta)
            pheromone_matrix = self.updatePheromone(pheromone_matrix, evaporation_rate, solutions)

            for solution in solutions:
                if solution.cost < best_cost:
                    best_cost = solution.cost
                    best_solution = solution
                    count += 1
                    generations_since_update = 0

        end_time = time.time()
        # Return the best solution found
        results['cost'] = best_solution.cost
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = best_solution
        return results

    def updatePheromone(self, pheromone_matrix, evaporation_rate, solutions):
        # Update pheromone matrix with evaporation and pheromone deposits
        pheromone_matrix = (1 - evaporation_rate) * pheromone_matrix

        for path in solutions:
            for i in range(len(path.route) - 1):
                pheromone_matrix[path.route[i]._index][path.route[i + 1]._index] += 1 / path.cost

        return pheromone_matrix


    def findAntsPaths(self, num_ants, pheromone_matrix, distance_matrix, alpha, beta):
        # Find the paths for each ant in the current iteration

        # Initialize empty lists for ants_paths and ants_costs
        solutions = []

        #  Helper function that will be executed by each thread
        def findPath(ant_index):
            # Find the path for the current ant
            solution = self.findAntPath(pheromone_matrix, distance_matrix, alpha, beta)
            solutions.append(solution)

        # Create and start a thread for each ant
        threads = []
        for i in range(num_ants):
            thread = threading.Thread(target=findPath, args=(i,))
            threads.append(thread)
            thread.start()

        # Wait threads to finish
        for thread in threads:
            thread.join()

        return solutions

    def findAntPath(self, pheromone_matrix, distance_matrix, alpha, beta):
        cities = self._scenario.getCities()
        current_node = random.randint(0, len(cities) - 1)  # Select a random node as the starting node
        path = []
        path_index = []

        # Uses set difference to remove the current node from the set of unvisited nodes
        unvisited_cities = set(range(len(distance_matrix))) - {current_node}
        path.append(cities[current_node])

        while unvisited_cities:
            # Select the next node based on the probability distribution
            next_node = self.calculateProbability(current_node, unvisited_cities, pheromone_matrix, distance_matrix,
                                                  alpha, beta)
            unvisited_cities.remove(next_node)
            path.append(cities[next_node])
            path_index.append(next_node)
            current_node = next_node
        solution = TSPSolution(path)
        return solution

    def calculateProbability(self, current_node, remaining_nodes, pheromone_matrix, distance_matrix, alpha, beta):
        if len(remaining_nodes) == 1:
            return next(iter(remaining_nodes))

        # Calculate the numerator for the probability formula by raising the pheromone values between the current node and each remaining node to the power of alpha.
        numerator = np.power(pheromone_matrix[current_node][list(remaining_nodes)], alpha)
        # Calculate the denominator for the probability formula by raising the distance values between the current node and each remaining node to the power of beta.
        denominator = np.power(distance_matrix[current_node][list(remaining_nodes)], beta)
        for i in range(len(denominator)):
            if denominator[i] == 0:
                denominator[i] = 10 ** -10

        probability = numerator / denominator
        probability = probability / np.sum(probability)  # Normalize the probability distribution so that it sums to 1
        for i in range(len(probability)):
            if np.isnan(probability[i]):
                probability[i] = 10 ** -10

        # Select the next node based on the probability distribution
        next_node = None
        try:
            next_node = np.random.choice(list(remaining_nodes), p=probability)  # find the value not the index
        except:
            probability = probability / np.sum(probability)
            next_node = np.random.choice(list(remaining_nodes), p=probability)
        return next_node
