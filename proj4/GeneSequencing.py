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

import random
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3
BAND_WIDTH = 7

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        self.MaxCharactersToAlign = None

        # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
        # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
        # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        SEQUENCE1 = "attggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccatgcccgtgggcttggtcttgtcatagtgctgacatttgtggttccttggtttttgttctctgccagtgacgtgtccattcggcgccagcagcccacccataggttgcataatggcaaagatgggcaaatacggtctcggcttcaaatgggccccagaatttccatggatgcttccgaacgcatcggagaagttgggtagccctgagaggtcagaggaggatgggttttgcccctctgctgcgcaagaaccaaaaactaaaggaaaaactttggttaatcacgtgagggtgaattgtagccggcttccagctttggaatgctgtgttcagtctgccataatccgtgatatttttgtagatgaggatccccagaaggtggaggcctcaactatgatggcattgcagttcggtagtgccgtcttggttaagccatccaagcgcttgtctattcaggcatggactaatttgggtgtgcttcccaaaacagctgccatggggttgttcaagcgcgtctgcctgtgtaacaccagggagtgctcttgtgacgcccacgtggcctttcacctttttacggtccaacccgatggtgtatgcctgggtaatggccgttttataggctggttcgttccagtcacagccataccggagtatgcgaagcagtggttgcaaccctggtccatccttcttcgtaagggtggtaacaaagggtctgtgacatccggccacttccgccgcgctgttaccatgcctgtgtatgactttaatgtagaggatgcttgtgaggaggttcatcttaacccgaagggtaagtactcctgcaaggcgtatgctcttcttaagggctatcgcggtgttaagcccatcctgtttgtggaccagtatggttgcgactatactggatgtctcgccaagggtcttgaggactatggcgatctcaccttgagtgagatgaaggagttgttccctgtgtggcgtgactccttggatagtgaagtccttgtggcttggcacgttgatcgagatcctcgggctgctatgcgtctgcagactcttgctactgtacgttgcattgattatgtgggccaaccgaccgaggatgtggtggatggagatgtggtagtgcgtgagcctgctcatcttctcgcagccaatgccattgttaaaagactcccccgtttggtggagactatgctgtatacggattcgtccgttacagaattctgttataaaaccaagctgtgtgaatgcggttttatcacgcagtttggctatgtggattgttgtggtgacacctgcgattttcgtgggtgggttgccggcaatatgatggatggctttccatgtccagggtgtaccaaaaattatatgccctgggaattggaggcccagtcatcaggtgttataccagaaggaggtgttctattcactcagagcactgatacagtgaatcgtgagtcctttaagctctacggtcatgctgttgtgccttttggttctgctgtgtattggagcccttgcccaggtatgtggcttccagtaatttggtcttctgttaagtcatactctggtttgacttatacaggagtagttg"
        SEQUENCE2 = "ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccatgcccgtgggcttggtcttgtcatagtgctgacatttgtggttccttggtttttgttctctgccagtgacgtgtccattcggcgccagcagcccacccataggttgcataatggcaaagatgggcaaatacggtctcggcttcaaatgggccccagaatttccatggatgcttccgaacgcatcggagaagttgggtagccctgagaggtcagaggaggatgggttttgcccctctgctgcgcaagaaccaaaaactaaaggaaaaactttgattaatcacgtgagggtggattgtagccggcttccagcattggagtgctgtgttcagtctgccataatccgtgatatttttgttgatgaggatcccttgaatgtggaggcctcaactatgatggcattgcagttcggtagtgctgtcttggtcaagccatccaagcgcttgtctattcaggcatgggctaagttgggtgtgctgcctaaaactccagccatggggttgttcaagcgcttctgcctgtgtaacaccagggagtgcgtttgtgacgcccacgtggcttttcaactttttacggtccagcctgatggtgtatgcctgggcaatggccgttttataggctggtttgtgccagtcacagccataccggcgtatgcgaagcagtggttgcaaccctggtccatccttcttcgtaagggtggtaacaaaggttctgtaacatccggccatttccgccgcgctgttaccatgcctgtgtatgactttaatgtggaggatgcttgtgaggaggttcatcttaacccgaagggtaagtactcccgcaaggcgtatgctcttcttaagggctatcgcggtgttaaatccatcctattcttggaccagtatggttgtgactatactgggcgtctcgccaagggtcttgaagactatggcgattgtactttggaagagatgaaggagttgtttcctgtgtggtgtgactccttggataatgaagttgttgtggcctggcatgttgatcgggatcctcgggctgttatgcgcttgcagactcttgctactgtacgttgcattggttatgtgggccaaccgaccgaggatttggttgatggagatgtggtagtgcgtgagcctgctcatctcctagcagccaatgccatcgtca"

        min_x = min(len(seq1), self.MaxCharactersToAlign)
        min_y = min(len(seq2), self.MaxCharactersToAlign)
        if banded is True:
            table = self.banded_algorithm(seq1, seq2)
            min_y = table[len(table) - 1][4]
            # min_y = self.get_last_item_index(table)
            score = min_y.get('cost')
        else:
            table = self.unrestricted_algorithm(SEQUENCE1, SEQUENCE2)
            score = table[min_x, min_y].get('cost')
        ###################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
        alignments = self.create_alignment(table, SEQUENCE1, SEQUENCE2)

        return {'align_cost': score, 'seqi_first100': alignments["seq1"], 'seqj_first100': alignments["seq2"]}

    def create_alignment(self, table, seq1, seq2):
        stack = []
        node = table[len(table) - 1, len(table[len(table) - 1]) - 1]
        stack.append(node)
        # Adds to the stack all the aligment
        while node["prev"] is not None:
            node = table[node["prev"][0], node["prev"][1]]
            stack.append(node)

        sequence1 = ""
        sequence2 = ""
        index1 = 0
        index2 = 0
        while len(stack) != 1 or index2 == 50 or index1 == 50:
            curNode = stack.pop()
            if stack[-1]['cost'] - curNode['cost'] == -3 or stack[-1]['cost'] - curNode['cost'] == 1:
                if index1 < len(seq1):
                    sequence1 += seq1[index1]
                if index2 < len(seq2):
                    sequence2 += seq2[index2]
                index1 += 1
                index2 += 1
            else:
                if stack[-1]['prev'][0] == curNode['prev'][0]:
                    if index1 < len(seq1):
                        sequence1 += seq1[index1]
                    sequence2 += '-'
                    index1 += 1
                else:
                    sequence1 += '-'
                    if index2 < len(seq2):
                        sequence2 += seq2[index2]
                    index2 += 1

        return {"seq1": sequence1, "seq2": sequence2}

    def unrestricted_algorithm(self, seq1, seq2):
        min_x = min(len(seq1), self.MaxCharactersToAlign)
        min_y = min(len(seq2), self.MaxCharactersToAlign)

        table = np.zeros((min_y + 1, min_x + 1), dtype=object)
        self.table_initializer(table)

        for i in range(1, len(seq2) + 1):
            for j in range(1, len(seq1) + 1):
                equalChar = seq1[j - 1] == seq2[i - 1]
                table[i, j] = self.min_distance(table[i, j - 1], table[i - 1, j - 1], table[i - 1, j], equalChar, i, j)

                # Check if align_length is exceeded
                if j == self.MaxCharactersToAlign:
                    break
            if i == self.MaxCharactersToAlign:
                break
        return table

    def banded_algorithm(self, seq1, seq2):
        min_y = min(len(seq1), self.MaxCharactersToAlign)
        # Create the first cases
        table = []
        self.banded_table_initializer(table)

        # Get bounds
        curRow = 1
        while curRow <= min_y:
            nextSpot = {"i": curRow, "j": 1}

            curColumn = max(1, curRow - MAXINDELS)
            upper_bound = min(curRow + MAXINDELS, len(seq2))

            while curColumn <= upper_bound:
                equalChar = seq2[curColumn - 1] == seq1[curRow - 1]
                # Append new row
                neighbours = self.calculateNeighbours(seq2, curRow, curColumn, table, nextSpot['j'])
                table[nextSpot['i']][nextSpot['j']] = self.banded_min_distance(neighbours['left'], neighbours['upLeft'],
                                                                               neighbours['up'], equalChar, nextSpot)
                curColumn += 1
                nextSpot['j'] += 1
            curRow += 1

            if min_y >= len(table) and curRow >= 4:
                table.append([0 for i in range(8)])
        return table

    def calculateNeighbours(self, seq2, curRow, curColumn, table, nextSpot):
        lower_bound = max(1, curRow - MAXINDELS)
        upper_bound = min(curRow + MAXINDELS, len(seq2))

        if upper_bound - lower_bound < 6 and curRow > 4:
            if curColumn == lower_bound:
                return {'left': 0, 'upLeft': table[curRow - 1][nextSpot], 'up': table[curRow - 1][nextSpot + 1]}
            elif curColumn == upper_bound:
                return {'left': table[curRow][nextSpot - 1], 'upLeft': table[curRow - 1][nextSpot],
                        'up': table[curRow - 1][nextSpot + 1]}
            else:
                return {'left': table[curRow][nextSpot - 1], 'upLeft': table[curRow - 1][nextSpot],
                        'up': table[curRow - 1][nextSpot + 1]}

        if curRow <= 4:
            return {'left': table[curRow][nextSpot - 1], 'upLeft': table[curRow - 1][nextSpot - 1],
                    'up': table[curRow - 1][nextSpot]}

        if curColumn == lower_bound:
            return {'left': 0, 'upLeft': table[curRow - 1][nextSpot], 'up': table[curRow - 1][nextSpot + 1]}
        elif curColumn == upper_bound:
            return {'left': table[curRow][nextSpot - 1], 'upLeft': table[curRow - 1][nextSpot], 'up': 0}
        else:
            return {'left': table[curRow][nextSpot - 1], 'upLeft': table[curRow - 1][nextSpot],
                    'up': table[curRow - 1][nextSpot + 1]}

    def min_distance(self, dist_left, dist_upLeft, dist_up, equalChar, i, j):
        if equalChar is False:
            min_dist = min(dist_left.get('cost'), dist_upLeft.get('cost'), dist_up.get('cost'))

            if min_dist == dist_left.get('cost'):
                min_location = dist_left
                min_dist += INDEL
            elif min_dist == dist_upLeft.get('cost'):
                min_location = dist_upLeft
                min_dist += SUB
            else:
                min_location = dist_up
                min_dist += INDEL
        else:
            min_location = dist_upLeft
            min_dist = dist_upLeft.get('cost') - 3

        # Tiebreakers
        if min_location == dist_left:
            if dist_upLeft.get('cost') == dist_left.get('cost') + SUB:
                min_location = dist_upLeft
                min_dist = dist_left.get('cost') + SUB
            elif dist_up.get('cost') == dist_left.get('cost') + INDEL:
                min_location = dist_up
                min_dist = dist_left.get('cost') + INDEL
        elif min_location == dist_up:
            if dist_upLeft.get('cost') == dist_up.get('cost') + SUB:
                min_location = dist_upLeft
                min_dist = dist_up.get('cost') + SUB

        min_location = min_location['loc']
        return {'prev': min_location, 'cost': min_dist, 'loc': (i, j)}

    def banded_min_distance(self, dist_left, dist_upLeft, dist_up, equalChar, nextSpot):
        if equalChar is False:
            # Check if within the band
            if dist_left != 0:
                left_dist = dist_left.get('cost')
            else:
                left_dist = float('inf')
            if dist_upLeft != 0:
                upleft_dist = dist_upLeft.get('cost')
            else:
                upleft_dist = float('inf')
            if dist_up != 0:
                up_dist = dist_up.get('cost')
            else:
                up_dist = float('inf')

            min_dist = min(left_dist + INDEL, upleft_dist + SUB, up_dist + INDEL)
            upLeftAdded = False
            if min_dist == upleft_dist + SUB:
                min_location = dist_upLeft
                upLeftAdded = True

            if min_dist == left_dist + INDEL and not upLeftAdded:
                min_location = dist_left

            if min_dist == up_dist + INDEL and not upLeftAdded:
                min_location = dist_up
        else:
            min_location = dist_upLeft
            min_dist = dist_upLeft.get('cost') - 3

        min_location = min_location['loc']
        return {'prev': min_location, 'cost': min_dist, 'loc': (nextSpot['i'], nextSpot['j'])}

    def table_initializer(self, table):
        # Populates the first column.
        curValue = 0
        for i in range(len(table)):
            table[i, 0] = {'prev': None, 'cost': curValue, 'loc': (i, 0)}
            curValue += INDEL

        curValue = 0
        # Populates the first row
        for j in range(len(table[0])):
            table[0, j] = {'prev': None, 'cost': curValue, 'loc': (0, j)}
            curValue += INDEL

    def banded_table_initializer(self, table):
        curValue = 0
        array = []
        for i in range(0, 8):
            if i < MAXINDELS + 1:
                array.append({'prev': None, 'cost': curValue, 'loc': (i, 0)})
                curValue += INDEL
            else:
                array.append(0)

        table.append(array)

        curValue = 5
        for j in range(1, MAXINDELS + 1):
            array = []
            for k in range(0, 8):
                if k == 0:
                    array.append({'prev': None, 'cost': curValue, 'loc': (0, j)})
                    curValue += INDEL
                else:
                    array.append(0)
            table.append(array)

    def get_last_item_index(self, table):
        last_item = None
        for i in range(len(table[len(table) - 1])):
            if table[len(table) - 1][i] == 0 and i != 0:
                return last_item
            else:
                last_item = table[len(table) - 1][i]

        return last_item
