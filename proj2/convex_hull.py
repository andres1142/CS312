from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25


#
# This is the class you have to complete.
#
class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        sortPoints = sorted(points, key=lambda point: point.x())
        # TODO: SORT THE POINTS BY INCREASING X-VALUE
        t2 = time.time()

        t3 = time.time()
        # this is a dummy polygon of the first 3 unsorted points
        polygon = self._divide_conquer(sortPoints, pause, view)
        # TODO: REPLACE THE LINE ABOVE WITH A CALL TO YOUR DIVIDE-AND-CONQUER CONVEX HULL SOLVER
        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        fullHull = [QLineF(polygon[i], polygon[(i + 1) % len(polygon)])
                    for i in range(len(polygon))]
        self.showHull(fullHull, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))

    def _divide_conquer(self, points, pause, view):
        numPoints = len(points)

        if numPoints == 1:
            return points
        leftHull = self._divide_conquer(points[:numPoints // 2], pause, view)
        rightHull = self._divide_conquer(points[numPoints // 2:], pause, view)

        # If there is only one point each just combine them
        if len(leftHull) == 1 and len(rightHull) == 1:
            leftHull.extend(rightHull)
            return leftHull

        # Find the right most point of the left hull and the left most point of the right hull
        # This is O(n)
        leftStart = leftHull.index(
            max(leftHull, key=lambda leftPoint: leftPoint.x()))
        rightStart = rightHull.index(
            min(rightHull, key=lambda rightPoint: rightPoint.x()))

        # Find the upper tangent
        # This is at worst O(n)
        i = leftStart
        j = rightStart
        left = True
        right = True
        slope = (rightHull[j].y() - leftHull[i].y()) / \
                (rightHull[j].x() - leftHull[i].x())
        while left or right:
            left = False
            right = False
            while True:
                newSlope = (rightHull[j].y() - leftHull[(i - 1) % len(leftHull)].y()) / (
                        rightHull[j].x() - leftHull[(i - 1) % len(leftHull)].x())
                if newSlope < slope:
                    left = True
                    slope = newSlope
                    i = (i - 1) % len(leftHull)
                else:
                    break
            while True:
                newSlope = (rightHull[(j + 1) % len(rightHull)].y() - leftHull[i].y()) / (
                        rightHull[(j + 1) % len(rightHull)].x() - leftHull[i].x())
                if newSlope > slope:
                    right = True
                    slope = newSlope
                    j = (j + 1) % len(rightHull)
                else:
                    break
        upper = (i, j)

        # Find the lower tangent
        # This is at worst O(n)
        i = leftStart
        j = rightStart
        left = True
        right = True
        slope = (rightHull[j].y() - leftHull[i].y()) / \
                (rightHull[j].x() - leftHull[i].x())
        while left or right:
            left = False
            right = False
            while True:
                newSlope = (rightHull[j].y() - leftHull[(i + 1) % len(leftHull)].y()) / (
                        rightHull[j].x() - leftHull[(i + 1) % len(leftHull)].x())
                if newSlope > slope:
                    left = True
                    slope = newSlope
                    i = (i + 1) % len(leftHull)
                else:
                    break
            while True:
                newSlope = (rightHull[(j - 1) % len(rightHull)].y() - leftHull[i].y()) / (
                        rightHull[(j - 1) % len(rightHull)].x() - leftHull[i].x())
                if newSlope < slope:
                    right = True
                    slope = newSlope
                    j = (j - 1) % len(rightHull)
                else:
                    break

        lower = (i, j)

        # Show recursion if selected
        if pause:
            self._show_recursion(leftHull, rightHull, upper, lower)

        # Combine the two hulls with upper and lower tangent
        # This is at worst O(n)
        final = []
        k = lower[0]
        final.append(leftHull[k])

        while k != upper[0]:
            k = (k + 1) % len(leftHull)
            final.append(leftHull[k])

        k = upper[1]
        final.append(rightHull[k])

        while k != lower[1]:
            k = (k + 1) % len(rightHull)
            final.append(rightHull[k])

        return final

    def _show_recursion(self, leftHull, rightHull, upper, lower):
        leftPrint = [QLineF(leftHull[i], leftHull[(i + 1) % len(leftHull)])
                     for i in range(len(leftHull))]
        rightPrint = [QLineF(rightHull[i], rightHull[(
                                                             i + 1) % len(rightHull)]) for i in
                      range(len(rightHull))]
        upperPrint = QLineF(leftHull[upper[0]], rightHull[upper[1]])
        lowerPrint = QLineF(leftHull[lower[0]], rightHull[lower[1]])
        self.showHull(leftPrint, RED)
        self.showHull(rightPrint, RED)
        self.showTangent([upperPrint, lowerPrint], BLUE)
        self.eraseHull(leftPrint)
        self.eraseHull(rightPrint)
        self.eraseTangent([upperPrint, lowerPrint])