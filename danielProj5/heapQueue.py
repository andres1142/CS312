from item import item
import math 

UNREACHABLE = -1
DEQUEUED = -2

class heapQueue:
    # initializes in O(n) time
    def __init__(self, length):
        self.queue = []
        self.keyMap = [UNREACHABLE for _ in range(length)]

    # O(1) time
    def isEmpty(self):
        return len(self.queue) == 0

    # O(log(n)) for the settleUp
    def insert(self, key, val):
        while len(self.keyMap) <= key:
            self.keyMap.append(-1)
        self.keyMap[key] = len(self.queue)
        self.queue.append(item(key, val))
        self.settleUp(self.keyMap[key])

    # O(log(n)) for the settleUp
    def decreaseKey(self, key, val):
        ind = self.keyMap[key]
        if ind == UNREACHABLE:
            self.insert(key, val)
            return
        self.queue[ind].val = val
        self.settleUp(ind)

    # O(log(n)) for the settleDown
    def deleteMin(self):
        minItem = self.queue[0]
        endInd = len(self.queue) - 1
        self.swap(0, endInd)
        self.keyMap[minItem.key] = DEQUEUED
        del self.queue[endInd]
        self.settleDown(0)
        return minItem.key

    # O(1) to swap
    def swap(self, ind1, ind2):
        item1 = self.queue[ind1]
        item2 = self.queue[ind2]
        self.queue[ind1] = item2
        self.queue[ind2] = item1
        self.keyMap[item1.key] = ind2
        self.keyMap[item2.key] = ind1

    # moves a node up to its correct place
    # at most log(n) swaps (from bottom to top of heap), and each swap takes O(1) time, so the total is O(log(n))
    def settleUp(self, ind):
        parent = self.parent(ind)
        if ind == 0 or self.queue[parent].val < self.queue[ind].val:
            return
        self.swapUp(ind)
        self.settleUp(parent)

    # moves a node down to its correct place
    # since heap is balanced (not deeper than log(n) layers), this comes out to O(log(n)) time
    def settleDown(self, ind):
        minInd = ind
        left = self.leftChild(ind)
        right = self.rightChild(ind)
        size = len(self.queue) - 1
        if left < size and self.queue[left].val < self.queue[minInd].val:
            minInd = left
        if right < size and self.queue[right].val < self.queue[minInd].val:
            minInd = right
        if minInd != ind:
            self.swap(ind, minInd)
            self.settleDown(minInd)

    ## HELPER FUNCTIONS ##
    # all O(1)
    def swapUp(self, ind):
        self.swap(ind, self.parent(ind))

    def parent(self, ind):
        return (ind - 1) // 2

    def leftChild(self, ind):
        return (2 * ind) + 1

    def rightChild(self, ind):
        return (2 * ind) + 2

    def size(self):
        return len(self.queue)
