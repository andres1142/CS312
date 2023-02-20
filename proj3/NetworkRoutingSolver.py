#!/usr/bin/python3


from CS312Graph import *
import time


class NetworkRoutingSolver:
    def __init__(self):
        self.network = None
        self.source = None
        self.dest = None
        self.distances = {}
        self.previous = {}

    def initializeNetwork(self, network):
        assert (type(network) == CS312Graph)
        self.network = network

    def getShortestPath(self, destIndex):
        assert(self.source is not None)

        self.dest = destIndex
        path_edges = []
        total_length = 0
        node = self.network.nodes[self.dest]
        while node.node_id != self.source:
            edge = self.previous[node.node_id]
            path_edges.append((edge.src.loc, edge.dest.loc, '{:.0f}'.format(edge.length)))
            total_length += edge.length
            node = edge.src
        path_edges.reverse()
        return {'cost': total_length, 'path': path_edges}

    def computeShortestPaths(self, srcIndex, use_heap=False):
        assert(self.network is not None)

        self.source = srcIndex
        self.distances = {node.node_id: float('inf') for node in self.network.nodes}
        self.distances[srcIndex] = 0
        self.previous = {node.node_id: None for node in self.network.nodes}
        visited = set()
        pq = PriorityQueue()

        t1 = time.time()

        pq.put(srcIndex, 0)

        while not pq.empty():
            current_node_id = pq.get()
            current_node = self.network.nodes[current_node_id]

            if current_node_id in visited:
                continue

            visited.add(current_node_id)

            for edge in current_node.neighbors:
                if edge.dest.node_id not in visited:
                    distance = self.distances[current_node_id] + edge.length
                    if distance < self.distances[edge.dest.node_id]:
                        self.distances[edge.dest.node_id] = distance
                        self.previous[edge.dest.node_id] = edge
                        pq.put(edge.dest.node_id, distance)

        t2 = time.time()
        return t2 - t1


class PriorityQueue:
    def __init__(self):
        self.queue = []

    def empty(self):
        return len(self.queue) == 0

    def put(self, item, priority):
        self.queue.append((item, priority))
        self.queue.sort(key=lambda x: x[1])

    def get(self):
        return self.queue.pop(0)[0]


class BinaryMinHeap:
    def __init__(self, network, dist_values):
        self.heap = [None,]
        self.end_position = 1
        self.root = 1
        self.indicies = {}
        for node in dist_values:
            self.insertNode(node, dist_values[node])

    def insertNode(self, node, length):
        if len(self.heap) == 1:
            self.heap.append((node, length))
            self.indicies[node] = 1
            return 1
        else:
            self.heap.append((node, length))
            self.bubbleUp((node, length), len(self.heap)-1)

    def getRightChild(self, position):
        if 2 * position + 1 < len(self.heap):
            return self.heap[2 * position + 1]
        return None

    def getLeftChild(self, position):
        if 2 * position < len(self.heap):
            return self.heap[2 * position]
        return None

    def getParent(self, position):
        return self.heap[position // 2]

    @staticmethod
    def getParentIndex(position):
        return position // 2

    @staticmethod
    def getRightChildIndex(position):
        return 2 * position + 1

    @staticmethod
    def getLeftChildIndex(position):
        return 2 * position

    def bubbleUp(self, x, i):
        p = i//2
        while i != 1 and self.heap[p][1] > x[1]:
            self.heap[i] = self.heap[p]
            self.indicies[self.heap[i][0]] = p
            self.indicies[self.heap[p][0]] = i
            i = p
            p = i // 2
        self.heap[i] = x
        self.indicies[self.heap[i][0]] = i

    def pop(self):
        return self.heap.pop()

    def decreaseKey(self, x, i):
        self.heap[self.indicies[x[0]]] = x
        self.bubbleUp(x, i)

    def deleteMin(self):
        if len(self.heap) == 1:
            return None
        else:
            min_value = self.heap[1]
            max_value = self.heap.pop()
            if len(self.heap) > 1:
                self.heap[1] = max_value
                self.siftDown(max_value, 1)
            return min_value

    def siftDown(self, x, i):
        c = self.minChild(i)
        while c != 0 and self.heap[c][1] < x[1]:
            self.heap[i] = (self.heap[c][0], self.heap[c][1])
            i = c
            c = self.minChild(i)
        self.heap[i] = x

    def minChild(self, index):
        if 2 * index > len(self.heap):
            return 0
        else:
            leftChild = self.getLeftChild(index)
            rightChild = self.getRightChild(index)
            if rightChild is None and leftChild is None:
                return 0
            if rightChild is None:
                return self.getLeftChildIndex(index)
            if leftChild[1] < rightChild[1]:
                return self.getLeftChildIndex(index)
            elif leftChild[1] > rightChild[1]:
                return self.getRightChildIndex(index)
            else:
                return self.getRightChildIndex(index)

    def __str__(self):
        return self.heap.__str__()
