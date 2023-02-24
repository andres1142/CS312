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
            if edge is None:
                return {'cost': float('inf'), 'path': path_edges}
            path_edges.append((edge.src.loc, edge.dest.loc, '{:.0f}'.format(edge.length)))
            total_length += edge.length
            node = edge.src
        path_edges.reverse()
        return {'cost': total_length, 'path': path_edges}

    def computeShortestPaths(self, srcIndex, use_heap=False):
        assert (self.network is not None)

        self.source = srcIndex
        self.distances = {node.node_id: float('inf') for node in self.network.nodes}
        self.distances[srcIndex] = 0
        self.previous = {node.node_id: None for node in self.network.nodes}
        visited = set()

        t1 = time.time()

        if use_heap:
            pq = BinaryMinHeap()
            pq.insertNode(srcIndex, 0)

            while pq.size() > 0:
                current_node_id, current_distance = pq.deleteMin()
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
                            pq.insertNode(edge.dest.node_id, distance)

        else:
            pq = PriorityQueue()
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
    def __init__(self):
        self.heap = []

    def insertNode(self, node, length):
        self.heap.append((node, length))
        self.bubbleUp(len(self.heap)-1)

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

    def bubbleUp(self, i):
        p = BinaryMinHeap.getParentIndex(i)
        while i != 0 and self.heap[p][1] > self.heap[i][1]:
            self.heap[i], self.heap[p] = self.heap[p], self.heap[i]
            i = p
            p = BinaryMinHeap.getParentIndex(i)

    def pop(self):
        if len(self.heap) == 0:
            return None
        elif len(self.heap) == 1:
            return self.heap.pop()
        else:
            min_value = self.heap[0]
            self.heap[0] = self.heap.pop()
            self.siftDown(0)
            return min_value

    def decreaseKey(self, x, i):
        self.heap[i] = x
        self.bubbleUp(i)

    def deleteMin(self):
        if len(self.heap) == 0:
            return None
        elif len(self.heap) == 1:
            return self.heap.pop()
        else:
            min_value = self.heap[0]
            max_value = self.heap.pop()
            self.heap[0] = max_value
            self.siftDown(0)
            return min_value

    def siftDown(self, i):
        c = self.minChild(i)
        while c != 0 and self.heap[c][1] < self.heap[i][1]:
            self.heap[i], self.heap[c] = self.heap[c], self.heap[i]
            i = c
            c = self.minChild(i)

    def minChild(self, index):
        left_child_index = BinaryMinHeap.getLeftChildIndex(index)
        right_child_index = BinaryMinHeap.getRightChildIndex(index)

        if right_child_index < len(self.heap):
            if self.heap[left_child_index][1] < self.heap[right_child_index][1]:
                return left_child_index
            else:
                return right_child_index
        elif left_child_index < len(self.heap):
            return left_child_index
        else:
            return 0

    def size(self):
        return len(self.heap)

    def __str__(self):
        return self.heap.__str__()