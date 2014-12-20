#!/usr/bin/env python

#############################################################################
##
## Copyright (C) 2014 Lukas Puchon
## All rights reserved.
##
## Program      : TransformFinder 3D
## Author       : Lukas Puchon
## Revision     : 1.0v
## Language     : SVK
## Platform     : Windows
##
##
##
## You may use this file under the terms of the BSD license as follows:
##
## Copyright (c) <year>, <copyright holder>
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of Author nor the names of its contributors
##       may be used to endorse or promote products derived from this software
##       without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
#############################################################################

#############################################################################
##
## Copyright (C) 2013 Riverbank Computing Limited.
## Copyright (C) 2010 Nokia Corporation and/or its subsidiary(-ies).
## All rights reserved.
##
## This file is part of the examples of PyQt.
##
## $QT_BEGIN_LICENSE:BSD$
## You may use this file under the terms of the BSD license as follows:
##
## "Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
##     the names of its contributors may be used to endorse or promote
##     products derived from this software without specific prior written
##     permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
## $QT_END_LICENSE$
##
#############################################################################

from __future__ import division
import os
import re

from numpy import *
from math import sqrt

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import (QPixmap, QIcon, QColor)

from ui_mainwindow import Ui_MainWindow


class MainWindow(QWidget):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        defaultFile1 = "proj.txt";
        defaultFile2 = "mer.txt";


        self.set1 = set()
        self.set2 = set()
        self.vstup1 = list()
        self.vstup2 = list()
        self.row2ID = list()
        self.presnost = 20

        self.currentDir = str(os.path.dirname(os.path.abspath(__file__))+"\\")

        self.allowedFiles = (".txt",".py", ".TXT")
        i=0
        self.allowedRegExpFiles = ""
        self.allowedOpenFiles = ""
        for fileType in self.allowedFiles:
            if i > 0:
                self.allowedRegExpFiles = self.allowedRegExpFiles+"|"
                self.allowedOpenFiles = self.allowedOpenFiles+" "
            self.allowedRegExpFiles = self.allowedRegExpFiles+"\\"+fileType
            self.allowedOpenFiles = self.allowedOpenFiles+"*"+fileType
            i=i+1
        
        
        self.ui = Ui_MainWindow()
	
        self.ui.setupUi(self)
        
        # comes in future release
        icoWindow = QPixmap("icon.png")
        icoCount = QPixmap("count.png")
        self.setWindowIcon(QIcon(icoWindow))
        self.ui.pushButtonPrepocitat.setIcon(QIcon(icoCount))
        self.ui.pushButtonUkazka.setVisible(False)
        self.ui.labelPresnost.setVisible(False)
        self.ui.spinBoxPresnost.setVisible(False)
        self.initTable()

        self.openAndReadPredloha(self.currentDir+defaultFile1)
        self.openAndReadMeranie(self.currentDir+defaultFile2)
        self.ui.lineEditSubor1.setText(self.currentDir+defaultFile1)
        self.ui.lineEditSubor2.setText(self.currentDir+defaultFile2)
        #self.on_pushButtonSubor1_released()
        #self.on_pushButtonSubor2_released()

        self.bod1 = False
        self.bod2 = False
        self.bod3 = False

    def on_tableWidgetZachytne_cellPressed(self, row, col):
        if self.ui.radioButton1bod.isChecked():
            self.ui.lineEdit1bod.setText(self.readRow(row))
            self.bod1 = True
        if self.ui.radioButton2bod.isChecked():
            self.ui.lineEdit2bod.setText(self.readRow(row))
            self.bod2 = True
        if self.ui.radioButton3bod.isChecked():
            self.ui.lineEdit3bod.setText(self.readRow(row))
            self.bod3 = True
        
    def read3points(self):
        bod1text = re.findall(r'[\w.-]+', self.ui.lineEdit1bod.text())
        bod2text = re.findall(r'[\w.-]+', self.ui.lineEdit2bod.text())
        bod3text = re.findall(r'[\w.-]+', self.ui.lineEdit3bod.text())
        return {'ids':(bod1text[0],bod2text[0],bod3text[0]), 
                'bod1':(float64(bod1text[1]),float64(bod1text[2]),float64(bod1text[3])),
                'bod2':(float64(bod2text[1]),float64(bod2text[2]),float64(bod2text[3])),
                'bod3':(float64(bod3text[1]),float64(bod3text[2]),float64(bod3text[3]))}
        
    def on_spinBoxPresnost_editingFinished(self):
        self.presnost = self.spinBoxPresnost.value()

    def on_lineEditSubor1_editingFinished(self):
        self.openAndReadPredloha(self.ui.lineEditSubor1.text())
        
    def on_lineEditSubor2_editingFinished(self):
        self.openAndReadMeranie(self.ui.lineEditSubor2.text())
        
    def on_pushButtonObnovit_released(self):
        self.initTable()
        
    def on_pushButtonVsetky_released(self):
        rows = self.ui.tableWidgetZachytne.rowCount()
        cols = self.ui.tableWidgetZachytne.columnCount()
        selection_model = self.ui.tableWidgetZachytne.selectionModel()
        zero = self.ui.tableWidgetZachytne.model().index(0, 0);
        bottomRight = self.ui.tableWidgetZachytne.model().index(rows - 1, cols - 1);
        selection = QItemSelection(zero, bottomRight)
        #selection.merge(QItemSelection(zero, bottomRight),QItemSelectionModel.Select)
        self.ui.tableWidgetZachytne.selectionModel().select(selection, QItemSelectionModel.Select)

    # open input files from dialog browsers
    def on_pushButtonSubor1_released(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Otvor Subor - Body Predlohy", self.currentDir, "Podporovane Subory ("+self.allowedOpenFiles+")")
        if fileName:
            self.openAndReadPredloha(fileName)
                
    def on_pushButtonSubor2_released(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Otvor Subor - Body z Merania", self.currentDir, "Podporovane Subory ("+self.allowedOpenFiles+")")
        if fileName:
            self.openAndReadMeranie(fileName)
                
    # drag and drop of input files
    
    def findFilePath(self, text):
        rx = QRegExp("((file:///).*("+self.allowedRegExpFiles+"))")
        rx.setPatternSyntax(QRegExp.RegExp2)
        index = rx.indexIn(text)
        if -1 != index:
            return (index, rx.cap(0)[8:])

        return (-1, "")
        
    def on_plainTextEditPredloha_textChanged(self):
        text = self.ui.plainTextEditPredloha.document().toPlainText()
        possibleFile = self.findFilePath(text)
        if possibleFile[0] != -1:
            self.openAndReadPredloha(possibleFile[1])
        else:
            di = self.parseTableTxt(text)
            self.set1 = di['set']
            self.initTable()
            
    def openAndReadPredloha(self, fileName):
        file = QFile(fileName)
        if not file.open(QFile.ReadOnly | QFile.Text):
            QMessageBox.warning(self, "Transformacia Bodov",
                    "Chyba pri otvarani suboru %s:\n%s." % (fileName, file.errorString()))
        else:
            di = self.parseTableTxt(self.readTextFile(file))
            self.ui.plainTextEditPredloha.setPlainText(di['text'])
            self.set1 = di['set']
            self.initTable()
            
    def parseTableTxt(self, text ):
        rows = text.split("\n")
        rows.sort()
        outText = ""
        setPts = set()
        prevRow = list()
        for row in rows:
            if len(row) > 6 and row != prevRow :
                prevRow = row
                outText = outText+row+"\n"
                cols = re.findall(r'[\w.-]+', row) # replaced row.split(" ") - working for multiple delimiters
                print(cols)
                if len(cols) > 3:
                    setPts.add((cols[0],cols[1],cols[2],cols[3]))

        return {'text':outText,'set':setPts}

    def readRow(self, row):
        tmpList = list()
        for col in range(0,self.ui.tableWidgetZachytne.columnCount()):
            item = self.ui.tableWidgetZachytne.item(row, col)
            for i in range(0,3):
                if not self.is_number(item.text()):
                    QMessageBox.warning(self, "Chyba vstupu", "Medzi zachytnymi bodmi v tabulke je neciselny vstup!\nRiadok: %s"  % item.row())
                    return False
                if i == item.column():
                    tmpList.append(float64(item.text()))
            
            if len(tmpList) == 3:
                return str(self.row2ID[row])+" "+str(tmpList[0])+" "+str(tmpList[1])+" "+str(tmpList[2])
                
    def initTable(self):
        commonIDs = set()
        set1IDs = set()
        set2IDs = set()
        for a in self.set1:
            set1IDs.add(a[0])

        for a in self.set2:
            set2IDs.add(a[0])

        commonIDs = set1IDs.intersection(set2IDs)
        commonRows = list()
        listSet2 = list(self.set2)
        for a in self.sorted_nicely( commonIDs ):
            for b in listSet2:
                if a == b[0]:
                    commonRows.append(b)

        # see dependent vector
        dependentVectors = dict()
        for a in commonRows:
            for b in commonRows:
                try:
                    if self.is_number(a[1]) and self.is_number(b[1]) and self.is_number(a[2]) and self.is_number(b[2]) and self.is_number(a[3]) and self.is_number(b[3]): 
                        kx = float(a[1]) / float(b[1])
                        ky = float(a[2]) / float(b[2])
                        kz = float(a[3]) / float(b[3])
                        if kx == ky and kx == kz and a != b:
                            randColor = QColor.fromHsv(qrand() % 256, 230, 220)
                            dependentVectors[a[0]]=randColor
                            dependentVectors[b[0]]=randColor
                except ZeroDivisionError:
                    xxx = 0
                    #print("division by zero")
        self.ui.tableWidgetZachytne.setRowCount(len(commonIDs))
        self.ui.tableWidgetZachytne.setColumnCount(3)
        rowIDs = list()
        colIDs = list(["X","Y","Z"])
        self.ui.tableWidgetZachytne.setHorizontalHeaderLabels(colIDs)
        self.ui.tableWidgetZachytne.setSelectionBehavior(QTableWidget.SelectRows)
        row = 0
        for curRow in commonRows:

            xItem = QTableWidgetItem(curRow[1])
            yItem = QTableWidgetItem(curRow[2])
            zItem = QTableWidgetItem(curRow[3])
            
            # if dependent : rows will be filled with same color
            if curRow[0] in dependentVectors:
                xItem.setBackground(dependentVectors[curRow[0]])
                yItem.setBackground(dependentVectors[curRow[0]])
                zItem.setBackground(dependentVectors[curRow[0]])
            # bg color with red if not a number
            if not self.is_number(curRow[1]):
                xItem.setBackground(Qt.red)
            # bg color with red if not a number
            if not self.is_number(curRow[2]):
                yItem.setBackground(Qt.red)
            # bg color with red if not a number
            if not self.is_number(curRow[3]):
                zItem.setBackground(Qt.red)

            self.ui.tableWidgetZachytne.setItem(row, 0, xItem)
            self.ui.tableWidgetZachytne.setItem(row, 1, yItem)
            self.ui.tableWidgetZachytne.setItem(row, 2, zItem)
            rowIDs.append(str(curRow[0]))
            row=row+1
            
        self.ui.tableWidgetZachytne.setVerticalHeaderLabels(rowIDs)
        self.row2ID = rowIDs
        
    def sorted_nicely(self, l ): 
        """ Sort the given iterable in the way that humans expect.""" 
        convert = lambda text: int(text) if text.isdigit() else text 
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        return sorted(l, key = alphanum_key)

    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    def on_plainTextEditMeranie_textChanged(self):
        text = self.ui.plainTextEditMeranie.document().toPlainText()
        possibleFile = self.findFilePath(text)
        if possibleFile[0] != -1:
            self.openAndReadMeranie(possibleFile[1])
        else:
            di = self.parseTableTxt(text)
            self.set2 = di['set']
            self.initTable()

    def openAndReadMeranie(self, fileName):
        file = QFile(fileName)
        if not file.open(QFile.ReadOnly | QFile.Text):
            QMessageBox.warning(self, "Transformacia Bodov",
                                "Chyba pri otvarani suboru %s:\n%s." % (fileName, file.errorString()))
        else:
            di = self.parseTableTxt(self.readTextFile(file))
            self.ui.plainTextEditMeranie.setPlainText(di['text'])
            self.set2 = di['set']
            self.initTable()
        
    def readTextFile(self, file):
        inf = QTextStream(file)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        text = inf.readAll()
        QApplication.restoreOverrideCursor()
        return text

    # saving of outputs
    def on_pushButtonUlozit_released(self):
        self.ui.plainTextEditTransformovane.appendPlainText("ulozit")
        if self.ui.plainTextEditTransformovane.document().isModified():
            self.saveAs()
            self.ui.plainTextEditTransformovane.document().setModified(False)
    
    def saveAs(self):
        fileName, _ = QFileDialog.getSaveFileName(self)
        if fileName:
            return self.saveFile(fileName)

        return False

    def saveFile(self, fileName):
        file = QFile(fileName)
        if not file.open(QFile.WriteOnly | QFile.Text):
            QMessageBox.warning(self, "Transformacia Bodov",
                    "Nemozno zapisat do suboru %s:\n%s." % (fileName, file.errorString()))
            return False

        outf = QTextStream(file)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        outf << self.ui.plainTextEditTransformovane.toPlainText()
        QApplication.restoreOverrideCursor()
        return True

    # future projects
    def on_pushButtonUkazka_released(self):
            QMessageBox.information(self, "Not Implemented Yet", "This should show all three datasets in opengl for visualization")
    
    # the big thing : MATH JOB
    def on_pushButtonPrepocitat_released(self):
        self.vstup1 = list()
        self.vstup2 = list()
        self.vstup2all = list()
        self.pts3id = list()
        tmpList = list()
        selectionIDs = list()
        counter = 0;
        for item in self.ui.tableWidgetZachytne.selectedItems():
            print("row %s" % item.row())
            for i in range(0,3):
                if not self.is_number(item.text()):
                    QMessageBox.warning(self, "Chyba vstupu", "Medzi zachytnymi bodmi v tabulke je neciselny vstup!\nRiadok: %s"  % item.row())
                    return False
                if i == item.column():
                    print(float64(item.text()))
                    tmpList.append(float64(item.text()))
            print("tmpList %s" % tmpList)
            if len(tmpList) == 3:
                self.vstup2.append(tmpList)
                selectionIDs.append(self.row2ID[item.row()])
                counter = counter + 1
                tmpList = list()
                for rowIn1 in self.set1:
                    if rowIn1[0] == self.row2ID[item.row()]:
                        self.vstup1.append( (float64(rowIn1[1]),float64(rowIn1[2]),float64(rowIn1[3])) )
                        break
        tmpList = list()
        for row in range(0,self.ui.tableWidgetZachytne.rowCount()):
            for col in range(0,self.ui.tableWidgetZachytne.columnCount()):
                item = self.ui.tableWidgetZachytne.item(row, col)
                for i in range(0,3):
                    if not self.is_number(item.text()):
                        QMessageBox.warning(self, "Chyba vstupu", "Medzi zachytnymi bodmi v tabulke je neciselny vstup!\nRiadok: %s"  % item.row())
                        return False
                    if i == item.column():
                        tmpList.append(float64(item.text()))
            
                if len(tmpList) == 3:
                    self.pts3id.append(self.row2ID[row])
                    self.vstup2all.append(tmpList)
                    tmpList = list()
                    
        pts1 = array(self.vstup1).reshape(counter,3)
        print("fixne body z predlohy")
        print(pts1)
        
        pts2 = array(self.vstup2).reshape(counter,3)
        print("fixne body z merania")
        print(pts2)

        pts2all = array(self.vstup2all).reshape(self.ui.tableWidgetZachytne.rowCount(),3)
        print("toto vsetko pretransformujeme")
        print(pts2all)
        
       # try:
        # be aware of transformation direction : from 1 to 2
        # in our case we want the other way so the points need to be swapped
        transformation = self.t_svd(pts2, pts1)
        pts3 = self.transformPoints(pts2all, transformation['rot'], transformation['trl'])
        print("vysledok")
        print(pts3)
        self.ui.lineEditChyba.setText(str(transformation['err']))
        cnt = 0
        for x in pts3:
            self.ui.plainTextEditTransformovane.appendPlainText(str(self.pts3id[cnt])+" "+str(x[0])+" "+str(x[1])+" "+str(x[2]))
            cnt = cnt +1

        if self.bod1 and self.bod2 and self.bod3:
            triBody = self.read3points()
            vstup1 = list()
            vstup2 = (triBody['bod1'],triBody['bod2'],triBody['bod3'])
            for bod in range(0,3):
                for rowIn1 in self.set1:
                    if rowIn1[0] == triBody['ids'][bod]:
                        vstup1.append( (float64(rowIn1[1]),float64(rowIn1[2]),float64(rowIn1[3])) )
                        break
            print("vstup1 pred alg 3 points :")
            print(vstup1)
            print("vstup2 pred alg 3 points :")
            print(vstup2)
            transform3points = self.alg_3_pts(vstup2, vstup1)
            pts3 = self.transformPoints2(pts2all, transform3points['rot1'],transform3points['rot2'], transform3points['trl'])
            self.ui.lineEditChyba_2.setText(str(transform3points['err']))
            cnt = 0
            for x in pts3:
                self.ui.plainTextEditTransformovane_2.appendPlainText(str(self.pts3id[cnt])+" "+str(x[0])+" "+str(x[1])+" "+str(x[2]))
                cnt = cnt +1

                
    def alg_3_pts(self, A, B):
        print ("A")
        print(A)
        print ("B")
        print(B)
        a1 = array(A[0])
        a2 = array(A[1])
        a3 = array(A[2])
        a_12 = a2 - a1
        print("a_12 = ")
        print(a_12)
        a_13 = a3 - a1
        print("a_13 = ")
        print(a_13)
        b1 = array(B[0])
        b2 = array(B[1])
        b3 = array(B[2])
        b_12 = b2 - b1
        b_13 = b3 - b1
        print("a1")
        print(a1)
        print("b1")
        print(b1)
        print ("translation, tl = b1 - a1")
        tl = b1 - a1
        print (tl)
        print("verify A + l, a1+l == b1")
        a1_ = a1 + tl
        a2_ = a2 + tl
        a3_ = a3 + tl
        print(a1_)
        print(a2_)
        print(a3_)
        print("\n----\nR =  rot_calc(a_12,b_12)")
        R = self.rot_calc(a2_,b2)
        R = self.rot_calc(a_12,b_12)
        print("verify R.(A-tl) , dot(R,a@) == b@")
        a1__ = dot(R,a1_)
        a2__ = dot(R,a2_)
        a3__ = dot(R,a3_)
        print(a1__)
        print(a2__)
        print(a3__)
        print("\n----\nR2 = rot_calc(dot(R,a_13),dot(R,b_13))")
        R2 = self.rot_calc(dot(R,a_13),b_13)
        print("verify R2+R.(A-tl) , dot(R,a@) == b@")
        a1___ = dot(R2,dot(R,a1_))
        a2___ = dot(R2,dot(R,a2_))
        a3___ = dot(R2,dot(R,a3_))
        print(a1___)
        print(a2___)
        print(a3___)
        mserr = 0.0
        dist = 0.0
        i =0
        for bod in (b1 - a1___, b2 - a2___, b3 - a3___):
            for i in range(0,3):
                print(i)
                x = bod[i]
                print(x)
                dist = dist + power(x,2)
                print(dist)
            
        mserr = sqrt(dist) / 3
        print("mserr === ")
        print(mserr)
        
        
        return {'trl':tl,'rot1':R,'rot2':R2, 'err':mserr}

    
    
    def rot_calc(self, a, b):
        print("vector a :")
        print(a)
        print("vector b :")
        print(b)
        print("unit vector au :")
        au = a/linalg.norm(a)
        print(au)
        print("unit vector bu :")
        bu = b/linalg.norm(b)
        print(bu)
        print("c = dot ( au, bu )")
        c = dot( au, bu )
        print (c)
        print ("v = cross ( a, b )")
        v = cross(au,bu,axis=0)
        print(v)
        print("s = ||v|| = norm ( v )")
        s = linalg.norm(v)
        print(s)
        print ("vx = skew symmetric cross product of vector v")
        vx = [[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]
        print(vx)
        print("R = (1-c)/power(s,2) * linalg.matrix_power(vx,2) + vx +  eye(3,3)")
        R = (1-c)/power(s,2) * linalg.matrix_power(vx,2) + vx +  eye(3,3)
        print(R)
        print("verify R*au == bu")
        print(dot(R,au))
        print("verify R*a =aprox= b")
        print(dot(R,a))
        return R

    def t_svd(self, A, B):

        # function from : http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/SVDalgorithm.ipynb
        # Ideally, three non-colinear markers placed on a moving rigid body is everything we need to describe
        # its movement (translation and rotation) in relation to a fixed coordinate system. However, in pratical
        # situations of human motion analysis, markers are placed on the soft tissue of a deformable body and this
        # generates artifacts caused by muscle contraction, skin deformation, marker wobbling, etc. In this situation,
        # the use of only three markers can produce unreliable results. It has been shown that four or more markers
        # on the segment followed by a mathematical procedure to calculate the 'best' rigid-body transformation taking
        # into account all these markers produces more robust results
        # (Söderkvist & Wedin 1993; Challis 1995; Cappozzo et al. 1997).
        
        Am = mean(A, axis=0)           # centroid of m1
        Bm = mean(B, axis=0)           # centroid of m2
        #print("Am")
        #print(Am)
        #print("Bm")
        #print(Bm)
        
        M = dot((B - Bm).T, (A - Am))  # considering only rotation
        #print("M")
        #print(M)
        
        # singular value decomposition - decomposes system of at least 3 lineary independent rows to orthonormal matrices
        U, S, Vt = linalg.svd(M)
        #print("U")
        #print(U)
        #print("Vt")
        #print(Vt)
        #print("S")
        #print(S)
        
        # rotation matrix
        R = dot(U, dot(diag([1, 1, linalg.det(dot(U, Vt))]), Vt))
        # translation vector
        L = B.mean(0)  - dot(R, A.mean(0))

        # RMSE from fixed points
        err = 0
        for i in range(A.shape[0]):
            Bp = dot(R, A[i, :]) + L
            err += sum((Bp - B[i, :])**2)
        RMSE = sqrt(err/A.shape[0]/3)
        #print("R")
        #print(R)
        #print ("L")
        #print (L)
        #print ("RMSE")
        #print (RMSE)
        
        return {'rot':R, 'trl':L, 'err':RMSE}
    
    def transformPoints(self, ptsIn, rotation, translation):
        count = ptsIn.shape[0]
        print("count")
        print(count)
        print("( ---------=  \n  vystupne bodiky")
        ptsOutList = list()
        for i in range(ptsIn.shape[0]):
            ptsOut = dot(rotation, ptsIn[i, :]) + translation
            ptsOutList.append(ptsOut)
            print (ptsOutList)
        return ptsOutList

    def transformPoints2(self, ptsIn, R, R2, translation):
        count = ptsIn.shape[0]
        print("count")
        print(count)
        print("( ---------=  \n  vystupne bodiky")
        ptsOutList = list()
        for i in range(ptsIn.shape[0]):
            ptsOut = dot(R2, dot(R, ptsIn[i, :] + translation)) 
            ptsOutList.append(ptsOut)
            print (ptsOutList)
        return ptsOutList
                
if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)
	
    mainwindow = MainWindow()
    mainwindow.show()
	
    sys.exit(app.exec_())
