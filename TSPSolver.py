#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4':
# 	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
from heapq import heappush, heappop


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None
        self.total = 0

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
        results['max'] = 0
        results['total'] = 0
        results['pruned'] = 0
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

    def branchAndBound(self, time_allowance=60.0):
        priority_queue = PriorityQueue()
        cities = self._scenario.getCities()
        ncities = len(cities)
        max_queue_size = 0
        # foundTour = False
        results = self.defaultRandomTour()
        print("Random tour bound: " + str(results['cost']))
        results['count'] = 0
        results['pruned'] = 0
        # at this point cost is equal to the random
        start_time = time.time()
        initial_state = self.generate_initial_state(cities)
        print("Initial bound: " + str(initial_state.bound))
        priority_queue.add(initial_state)
        self.total = 1
        print(initial_state.matrix)
        while time.time() - start_time < time_allowance and priority_queue.get_size() > 0:
            if priority_queue.get_size() > max_queue_size:
                max_queue_size = priority_queue.get_size()
            state = priority_queue.pop_off()

            # check if this state is a complete path
            if len(state.path) == ncities:
                if state.bound < results['cost']:
                    cities_path = [cities[state.path[j]] for j in range(len(state.path))]
                    results['soln'] = TSPSolution(cities_path)
                    results['bssf'] = TSPSolution(cities_path)
                    results['cost'] = state.bound
                    results['count'] += 1
                    print("Found a solution")
            # check if is bssf
            elif state.bound >= results['cost']:
                results['pruned'] += 1
                print("Pruning...")
            else:
                # child and everything
                children = self.expand_tree(state)

                # if it is, then find children.
                # if child bound < bssf, put on priority queue
                # if not, prune.

                for i in range(len(children)):
                    if children[i].bound < results['cost']:
                        priority_queue.add(children[i])
                        print("Adding to priority queue")
                    else:
                        results['pruned'] += 1
                        print("Pruning...")
        results['cost'] = results['soln'].cost
        end_time = time.time()
        results['time'] = end_time - start_time
        results['max'] = max_queue_size
        results['total'] = self.total
        return results

    def expand_tree(self, state):
        new_states = []
        for i in range(0, len(state.matrix)):
            if i not in state.path:
                new_bound = state.bound + state.matrix[state.path[-1]][i]
                if new_bound < math.inf:
                    new_matrix = state.matrix.copy()
                    #set row and column to infinity
                    for j in range(len(new_matrix)):
                        new_matrix[state.path[-1]][j] = math.inf
                        new_matrix[j][i] = math.inf
                    new_matrix[i][state.path[-1]] = math.inf

                    #Create a new state from the old matrix
                    new_path = state.path.copy()
                    new_path.append(i)
                    new_state = State(new_matrix, new_path, new_bound)

                    #calculate bound and add to list
                    new_state = self.calculate_bound(new_state)
                    new_states.append(new_state)
                    self.total += 1

        return new_states

    def generate_initial_state(self, cities):
        #create matrix
        matrix = np.zeros((len(cities), len(cities)))
        #initialize costs
        for i in range(len(cities)):
            for j in range(len(cities)):
                if i == j:
                    matrix[i][j] = math.inf
                matrix[i][j] = cities[i].costTo(cities[j])
        path = [0]
        bound = 0
        # create state and calculate the bound
        state = State(matrix, path, bound)
        state = self.calculate_bound(state)
        return state

    def calculate_bound(self, state):
        bound = state.bound
        # matrix = state.matrix
        n_rows = len(state.matrix)
        n_cols = len(state.matrix[0])

        # subtract the minimum value in each row from all cells in that row
        for i in range(n_rows):
            row_min = min(state.matrix[i])
            if row_min > 0 and row_min < math.inf:
                for j in range(n_cols):
                    state.matrix[i][j] -= row_min
                bound += row_min

        # subtract the minimum value in each column from all cells in that column
        for j in range(n_cols):
            col_min = min(state.matrix[i][j] for i in range(n_rows))
            if col_min > 0 and col_min < math.inf:
                for i in range(n_rows):
                    state.matrix[i][j] -= col_min
                bound += col_min
        state.bound = bound
        return state

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


class State:
    def __init__(self, matrix, path, bound):
        self.matrix = matrix
        self.path = path
        self.bound = bound


class PriorityQueue:
    def __init__(self):
        self.queue = []

    def add(self, state):

        depth = len(state.path) - 1

        # if there is no list at that depth, create a new one
        # with that state and add it to that depth
        if depth > len(self.queue) - 1:
            nodes = [state]
            self.queue.insert(depth, nodes)
        else:
            index = 0
            # while index is less than the length of the list at that depth
            while index < len(self.queue[depth]) and state.bound > self.queue[depth][index].bound:
                index += 1

            self.queue[depth].insert(index, state)

    def pop_off(self):
        # pop off at the deepest level with the lowest bound
        depth = len(self.queue) - 1
        for i in range(len(self.queue)):
            if len(self.queue[depth - i]) != 0:
                state = self.queue[depth - i].pop(0)
                return state

    def get_size(self):
        size = 0
        for i in range(len(self.queue)):
            size += len(self.queue[i])
        return size
