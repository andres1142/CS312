#!/usr/bin/python3
from math import inf

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

        min_x = min(len(seq2), self.MaxCharactersToAlign)
        min_y = min(len(seq1), self.MaxCharactersToAlign)

        if banded is True:
            table = self.banded_algorithm(seq1, seq2)
            min_y = (table[len(table) - 1][4], 4) if min_x > 100 else self.get_last_item_index(table)
            score = min_y
            alignments = self.create_banded_alignment(table, seq1, seq2, min_y)
        else:
            table = self.unrestricted_algorithm(seq1, seq2)
            score = table[min_x, min_y]
            alignments = self.create_alignment(table, seq1, seq2)

        sequence1 = alignments["seq1"][:100] if len(seq1) > 100 else alignments["seq1"]
        sequence2 = alignments["seq2"][:100] if len(seq2) > 100 else alignments["seq2"]

        return {'align_cost': score if banded is False else score[0], 'seqi_first100': sequence1,
                'seqj_first100': sequence2}

    def create_banded_alignment(self, table, seq1, seq2, lastColumn):

        if len(seq1) < self.MaxCharactersToAlign or len(seq2) < self.MaxCharactersToAlign:
            return {"seq1": "No Alignment Possible.", "seq2": "No Alignment Possible."}
        row = table.shape[0] - 1
        col = lastColumn[1]

        sequence1 = " "
        sequence2 = " "
        node = table[row][col]
        index = self.MaxCharactersToAlign
        while row != 0:
            equalChar = seq2[row - 1] == seq1[index]

            if row > 4:
                if col != 0 and col != 7:
                    up_cell = (table[row - 1, col + 1], [row - 1, col + 1])
                    left_cell = (table[row, col - 1], [row, col - 1])
                    upLeft_cell = (table[row - 1, col], [row - 1, col])

                # Work with lower bound
                elif col == 0:
                    up_cell = (table[row - 1, col + 1], [row - 1, col + 1])
                    left_cell = (float(inf), float(inf))
                    upLeft_cell = (table[row - 1, col], [row - 1, col])
                # Work with upper bound
                elif col == 7:
                    up_cell = (float(inf), float(inf))
                    left_cell = (table[row, col - 1], [row, col - 1])
                    upLeft_cell = (table[row - 1, col], [row - 1, col])

                min_value = min(up_cell[0] + INDEL, left_cell[0] + INDEL,
                                (upLeft_cell[0] + MATCH) if equalChar is True else (
                                            upLeft_cell[0] + SUB))


                if min_value == left_cell[0] + INDEL:
                    sequence1 = seq1[index - 1] + sequence1
                    sequence2 = "-" + sequence2
                    index -= 1
                    col = left_cell[1][1]
                    continue

                if min_value == up_cell[0] + INDEL:
                    sequence1 = "-" + sequence1
                    sequence2 = seq2[row - 1] + sequence2
                    row -= 1
                    col = up_cell[1][1]
                    continue

                if min_value == upLeft_cell[0] + MATCH or min_value == upLeft_cell[0] + SUB:
                    sequence1 = seq1[index - 1] + sequence1
                    sequence2 = seq2[row - 1] + sequence2
                    col = upLeft_cell[1][1]
                    row -= 1
                    index -= 1
                    continue
            else:

                up_cell = (table[row - 1,col], [row - 1, col])
                left_cell = (table[row, col - 1], [row, col - 1])
                upLeft_cell = (table[row - 1, col - 1], [row - 1, col - 1])

                min_value = min(up_cell[0] + INDEL, left_cell[0] + INDEL,
                                (upLeft_cell[0] + MATCH) if equalChar is True else (
                                        upLeft_cell[0] + SUB))

                if min_value == left_cell[0] + INDEL:
                    sequence1 = seq1[index - 1] + sequence1
                    sequence2 = "-" + sequence2
                    index -= 1
                    col = left_cell[1][1]
                    continue

                if min_value == up_cell[0] + INDEL:
                    sequence1 = "-" + sequence1
                    sequence2 = seq2[row - 1] + sequence2
                    row -= 1
                    col = up_cell[1][1]
                    continue

                if min_value == upLeft_cell[0] + MATCH or min_value == upLeft_cell[0] + SUB:
                    sequence1 = seq1[index - 1] + sequence1
                    sequence2 = seq2[row - 1] + sequence2
                    col = upLeft_cell[1][1]
                    row -= 1
                    index -= 1
                    continue



        return {"seq1": sequence1, "seq2": sequence2}

    def create_alignment(self, table, seq1, seq2):
        row = table.shape[0] - 1
        col = table.shape[1] - 1

        sequence1 = ""
        sequence2 = ""

        while row != 0 and col != 0:
            equalChar = seq2[row - 1] == seq1[col - 1]
            up_cell = table[row - 1, col]
            left_cell = table[row, col - 1]
            upLeft_cell = table[row - 1, col - 1]
            min_value = min(up_cell + INDEL, left_cell + INDEL,
                            (table[row - 1, col - 1] + MATCH) if equalChar is True else (table[row - 1, col - 1] + SUB))

            if min_value == left_cell + INDEL:
                sequence1 = seq1[col - 1] + sequence1
                sequence2 = "-" + sequence2
                col -= 1
                continue

            if min_value == up_cell + INDEL:
                sequence1 = "-" + sequence1
                sequence2 = seq2[row - 1] + sequence2
                row -= 1
                continue

            if min_value == upLeft_cell + MATCH or min_value == upLeft_cell + SUB:
                sequence1 = seq1[col - 1] + sequence1
                sequence2 = seq2[row - 1] + sequence2
                row -= 1
                col -= 1
                continue

        return {"seq1": sequence1, "seq2": sequence2}




    def unrestricted_algorithm(self, seq1, seq2):
        min_y = min(len(seq1), self.MaxCharactersToAlign)
        min_x = min(len(seq2), self.MaxCharactersToAlign)

        table = np.zeros((min_x + 1, min_y + 1), dtype=int)
        self.table_initializer(table)

        for i in range(1, min_y + 1):
            for j in range(1, min_x + 1):
                equalChar = seq2[j - 1] == seq1[i - 1]
                table[j, i] = self.min_distance(table[j, i - 1], table[j - 1, i - 1], table[j - 1, i], equalChar)
        return table

    def banded_algorithm(self, seq1, seq2):
        min_y = min(len(seq1), self.MaxCharactersToAlign)
        # Create the first cases
        table = self.banded_table_initializer()

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
                                                                               neighbours['up'], equalChar)
                curColumn += 1
                nextSpot['j'] += 1
            curRow += 1

            if min_y >= len(table) and curRow >= 4:
                table = np.vstack([table, np.zeros(table.shape[1], dtype=int)])
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


    def banded_min_distance(self, dist_left, dist_upLeft, dist_up, equalChar):
        if equalChar is False:
            # Check if within the band
            if dist_left != 0:
                left_dist = dist_left
            else:
                left_dist = 0
            if dist_upLeft != 0:
                upleft_dist = dist_upLeft
            else:
                upleft_dist = 0
            if dist_up != 0:
                up_dist = dist_up
            else:
                up_dist = 0

            min_dist = min(left_dist + INDEL, upleft_dist + SUB, up_dist + INDEL)
        else:
            min_dist = dist_upLeft + MATCH

        return min_dist

    def min_distance(self, dist_left, dist_upLeft, dist_up, equalChar):
        if equalChar is False:
            return min(dist_left + INDEL, dist_upLeft + SUB, dist_up + INDEL)
        else:
            return min(dist_left + INDEL, dist_upLeft + MATCH, dist_up + INDEL)

    def table_initializer(self, table):
        # Populates the first column.
        curValue = 0
        for i in range(len(table)):
            table[i, 0] = curValue
            curValue += INDEL

        curValue = 0
        # Populates the first row
        for j in range(len(table[0])):
            table[0, j] = curValue
            curValue += INDEL

    def banded_table_initializer(self):
        table = np.empty((0, 8), dtype=int)
        curValue = 0
        array = np.zeros(8, dtype=int)

        for i in range(0, 8):
            if i < MAXINDELS + 1:
                array[i] = curValue
                curValue += INDEL

        table = np.vstack([table, array])

        curValue = 5
        for i in range(3):
            array = np.zeros(8, dtype=int)
            array[0] = curValue
            table = np.vstack([table, array])
            curValue += 5
        return table

    def get_last_item_index(self, table):
        last_item = None
        for i in range(len(table[len(table) - 1])):
            if table[len(table) - 1][i] == 0 and i != 0:
                return (last_item, i - 1)
            else:
                last_item = table[len(table) - 1][i]
