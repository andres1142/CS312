#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

DIAGONAL = 0
LEFT = 1
UP = -1


# class Matrix:
#
# 	def __init__(self):
#


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment
    def printList(self, results, sequence1, sequence2, seq1Length, seq2Length):
        adjustedLength1 = seq1Length + 1
        adjustedLength2 = seq2Length + 1
        ##print('adj len 1:',adjustedLength1,' ad len 2:' , adjustedLength2)
        # print out the letters of the first sequence along the top
        print(end='- - ')
        for u in range(seq1Length):
            print(sequence1[u], end=' ')
        print('')
        # for each row
        for w in range(seq2Length + 1):
            # first row doesn't match a letter
            if w == 0:
                print('-', end=' ')
            else:
                # print out the associated letter at the front of the row
                print(sequence2[w - 1], end=' ');
            for x in range(seq1Length + 1):
                index = (w * adjustedLength1) + x
                print(results[index], end=' ')
            # print(results[index], end = ' ')
            print('')
        print('')
        print('')

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        toReturn = []

        # for each pair of sequences
        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                if (j <= i):
                    print("looking at sequence ", i, " and sequence ", j)
                    if banded:
                        hi = 0
                    else:
                        sequence1 = sequences[i]
                        sequence2 = sequences[j]
                        seq1Length = len(sequences[i])
                        seq2Length = len(sequences[j])
                        adjustedLength1 = seq1Length + 1
                        adjustedLength2 = seq2Length + 1
                        # print('adjustedLength1 =' ,adjustedLength1,' adjusted length2 = ',adjustedLength2)
                        # sequence 1 is the columns/spans the top of the table
                        # sequence 2 is the row / spans the side of the table
                        # make all cells have a value 0
                        # results = [[0 for x in range(seq1Length+1)] for y in range(seq2Length+1)]
                        # pointers = [[0 for x in range(seq1Length+1)] for y in range(seq2Length+1)]

                        # columns = [None]*(seq1Length+1)
                        # results = [columns]*(seq2Length+1)

                        # although I am using only one list, I treat is as a matrix.
                        # after having gone the length of seq, the next 'row' begins
                        results = (seq2Length + 1) * (seq1Length + 1) * [None]
                        # print('results len = ', len(results))
                        cost = 0
                        # populate 1st row with multiples of five
                        print('populating 1st row')
                        for k in range(seq1Length + 1):
                            results[k] = cost
                            cost += 5

                        cost = 5
                        # populate 1st column with multiples of 5
                        # start at the beginning of the second row
                        # print('populating 1st column')
                        index = seq1Length + 1
                        for k in range(seq2Length):
                            # print('line 96 index = ', index)
                            results[index] = cost
                            cost += 5
                            # move to the next "row"
                            index += (seq1Length + 1)

                        # calculate the cost for the rest
                        for l in range(1, adjustedLength2):  # 'rows'
                            for m in range(1, adjustedLength1):  # 'columns'
                                index = (l * adjustedLength1) + m
                                # print('index = ', index)
                                left = (results[index - 1]) + 5
                                up = (results[index - adjustedLength1]) + 5
                                diagonal = 200
                                # if the letters match
                                previousdiagIndex = index - adjustedLength1 - 1
                                if (sequence2[l - 1] == sequence1[m - 1]):
                                    diagonal = results[previousdiagIndex] - 3
                                else:
                                    diagonal = results[previousdiagIndex] + 1
                                # find the min
                                # diagonal preceeds up preceeds left
                                if diagonal < left and diagonal < up:
                                    results[index] = diagonal
                                    pointers[index] = DIAGONAL
                                elif up < diagonal and up < left:
                                    results[index] = up
                                    pointers[index] = UP
                                else:
                                    results[index] = left
                                    pointers[index] = LEFT

                        if i < 2 and j < 2:
                            self.printList(results, sequence1, sequence2, seq1Length, seq2Length)

                    score = results[((seq2Length + 1) * (seq1Length + 1)) - 1]
                    table.item(i, j).setText('{}'.format(int(score)))
                    table.repaint()
                    index = ((seq2Length + 1) * (seq1Length + 1)) - 1
                    # backtrace
                    aligment1 = []
                    aligment2 = []
                    seq1Index = seq1Length - 1
                    seq2Index = seq2Length - 1
                    while (index != 0):
                        temp = pointers[index]
                        if (seq1Index < 0):
                            seq1Index = 0
                        if (temp == DIAGONAL):
                            aligment1.append(sequence1[seq1Index])
                            seq1Index = seq1Index - 1
                            aligment2.append(sequence2[seq2Index])
                            seq2Index -= 1
                            index -= (adjustedLength1 + 1)

                        elif (temp == UP):
                            aligment2.append(sequence2[seq2Index])
                            aligment1.append('-')
                            seq2Index -= 1
                            index -= adjustedLength1
                        else:
                            aligment1.append(sequence1[seq1Index])
                            aligment2.append('-')
                            seq1Index = seq1Index - 1
                            index -= 1
                    string1 = ''.join(reversed(aligment1))
                    string2 = ''.join(reversed(aligment2))
                    if (j == 9 and i == 2):
                        print(string1)
                        print(string2)
                    s = {'align_cost': score, 'seqi_first100': string1, 'seqj_first100': string2}
                    jresults.append(s)
                    print('length of jresults = ', len(jresults))
                toReturn.append(jresults)
        return toReturn

# else:
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 					score = i+j;
# 					alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
# 						len(sequences[i]), align_length, ',BANDED' if banded else '')
# 					alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
# 						len(sequences[j]), align_length, ',BANDED' if banded else '')
# ###################################################################################################
# 					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
# 					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
# 					table.repaint()
# jresults.append(s)
# results.append(jresults)
# return results


# -----------------------Deprecated using two lists ----------------------------#
# for k in range(seq1Length+1):
# 	results[0][k] = cost
# 	cost += 5

# cost = 0
# #populate 1st column with multiples of 5
# print('populating 1st column')
# for k in range(seq2Length+1):
# 	results[k][0] = cost
# 	cost +=5


# calculate cost for the rest
# print('populating the rest')
# for l in range(1,seq2Length+1):
# 	for m in range(1,seq1Length+1):
# 		left = (results[l][m-1])+5
# 		up = (results[l-1][m])+5
# 		diagonal = 0
# 		#is it a match of sub?
# 		if(sequence2[l-1] == sequence1[m-1]):
# 			diagonal = results[l-1][m-1] - 3
# 		else:
# 			diagonal = results[l-1][m-1] + 1
# 		if left < up and left < diagonal:
# 			results[l][m] = left
# 		elif up  < left and up < diagonal:
# 			results[l][m] = up
# 		else:
# 			results[l][m] = diagonal


# previousRow = adjustedLength1 * [None]			# save the previous row for quicker lookup
# currentRow = adjustedLength1 * [None]					# save this for quicker lookup
# previousValue = math.inf
# #print('adjustedLength1 = ', adjustedLength1)
# for l in range(0,adjustedLength2):			#  'rows'
# 	for m in range(0,adjustedLength1):		#  'columns'
# 		index = (l*adjustedLength1)+m
# 		if(l == 0):
# 			#populate first row
# 			results[index] = cost
# 			currentRow[m] = cost
# 			cost+=5
# 		else:
# 			if(l == 1):
# 				cost = 5
# 			#populate first column
# 			if(m == 0):
# 				results[index] = cost
# 				currentRow[m] = cost
# 				previousValue = math.inf
# 				cost+=5
# 			else:
# 				print('current row = ', l, ' and previous row starts with ', previousRow[0])
# 				left = previousValue +5
# 				up = previousRow[m] + 5
# 				diagonal = previousRow[m-1]
# 				if(sequence2[l-1] == sequence1[m-1]):
# 					diagonal = diagonal-3
# 				else:
# 					diagonal = diagonal +1
# 				if diagonal < left and diagonal < up:
# 					currentRow[m] = diagonal
# 					previousValue = diagonal
# 					results[index] = diagonal
# 				elif up < diagonal and up < left:
# 					currentRow[m] = up
# 					previousValue = up
# 					results[index] = up
# 				else:
# 					results[index] = left
# 					currentRow[m] = left
# 					previousValue = left
# 	print('changing current row to previous row')
# 	previousRow = currentRow								# finished a row, previous row becomes current
# 	print('previous row index 0 is now:', previousRow[0])
# 	hey = [0,1,2]
# 	current = hey
# 	hey[0] = 10
# 	print('hlkjhlkjhl   ',current[0])
