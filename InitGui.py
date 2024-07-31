# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2024 - Harry van Langen <hvlanalysis@gmail.com>        *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

__Title__ = "fcSC"
__Author__ = "HarryvL"
__Url__ = "https://github.com/HarryvL/fcVM-workbench"
__Version__ = "1.1.0"
__Date__ = "2024/06/25"
__Comment__ = "first release"
__Forum__ = "https://forum.freecad.org/viewtopic.php?t=85474"
__Status__ = "initial development"
__Requires__ = "freecad version 0.19 or higher"
__Communication__ = "https://forum.freecad.org/viewtopic.php?t=85474"

import os
import sys
import dummySC
import FreeCAD
import FreeCADGui
from PySide2 import QtWidgets, QtGui, QtCore
from PySide2.QtWidgets import QTableWidgetItem

global FCmw
FCmw = FreeCADGui.getMainWindow()
global dir_name
dir_name = os.path.dirname(dummySC.file_path())
global double_validator
double_validator = QtGui.QDoubleValidator()
global int_validator
int_validator = QtGui.QIntValidator()
global res_show
res_show = FreeCAD.getHomePath() + "Mod/Fem/Resources/ui/ResultShow.ui"
global source_code_path
source_code_path = os.path.join(dir_name, 'source code')

sys.path.append(source_code_path)


class fcSCWorkbench(Workbench):
    Icon = os.path.join(dir_name, "icons", "fcFEM.svg")
    MenuText = "fcSC workbench"
    ToolTip = "Simplified plastic collapse analysis of Structural Concrete"

    def __init__(self):
        import task_result_mechanical
        sys.modules["femtaskpanels.task_result_mechanical"] = sys.modules[task_result_mechanical.__name__]

    def GetClassName(self):
        return "Gui::PythonWorkbench"

    def Initialize(self):
        from PySide2 import QtGui
        self.appendToolbar("fcSC", [])
        self.appendMenu("fcSC", [])
        self.palette_warning = QtGui.QPalette()
        self.palette_warning.setColor(QtGui.QPalette.Base, QtGui.QColor("orange"))
        self.palette_standard = QtGui.QPalette()
        self.palette_standard.setColor(QtGui.QPalette.Base, QtGui.QColor("white"))
        # self.palette.setColor(QtGui.QPalette.Text, QtGui.QColor("red"))

    def Activated(self):
        from PySide2 import QtCore
        global fcSC_window

        import dummySC
        self.dir_name = os.path.dirname(dummySC.file_path())

        self.doc = FreeCAD.activeDocument()
        if self.doc == None:
            self.doc = FreeCAD.newDocument("fcSC")
        self.file_name = self.doc.Label
        self.macro_file_path = os.path.join(self.dir_name, "source code", "fcSC.FCMacro")
        self.sum_file_path = os.path.join(self.dir_name, "source code", "fcSC_sum.FCMacro")
        self.disp_option = "incremental"
        self.eps_option = "eps_c"
        self.averaged_option = "unaveraged"
        self.model_option = "1"
        self.mat_obj = {}
        self.objID = {}
        self.deleted_objs = []

        class DocObserver(object):  # document Observer
            def __init__(self, workbench_instance):
                self.workbench_instance = workbench_instance

            def slotDeletedObject(self, obj):
                #
                # mat_obj = {App::MaterialObjectPython : [rhox, rhoy, rhoz, dx, dy, dz]}
                #
                # objID = {row : App::MaterialObjectPython}
                #
                self.workbench_instance.deleted_objs.append(obj.Label)
                for row in self.workbench_instance.objID:
                    if self.workbench_instance.objID[row] == obj:
                        del self.workbench_instance.mat_obj[obj]
                self.workbench_instance.check_mat_objs()
                self.workbench_instance.update_material()

            def slotChangedObject(self, obj, prop):
                self.workbench_instance.update_material()

            def slotActivateDocument(self, doc):
                if FreeCAD.activeDocument().Label[0:7] != "Unnamed":
                    self.workbench_instance.save_clicked()
                    self.workbench_instance.file_name = FreeCAD.activeDocument().Label
                    self.workbench_instance.mat_obj = {}
                    self.workbench_instance.objID = {}
                    self.workbench_instance.deleted_objs = []
                    self.workbench_instance.check_mat_objs()
                    self.workbench_instance.open_file()
                    self.workbench_instance.update_material()

                    fcSC_window.eps_c.setText("0.000")
                    fcSC_window.eps_s.setText("0.000")

            def slotFinishSaveDocument(self, doc, prop):
                self.workbench_instance.save_clicked()  # save under old file name
                self.workbench_instance.file_name = doc.Label
                self.workbench_instance.save_clicked()  # save under new file name

        self.obs = DocObserver(self)
        FreeCAD.addDocumentObserver(self.obs)

        ui_Path = os.path.join(self.dir_name, "user_interface", "fcSC.ui")

        fcSC_window = FreeCADGui.PySideUic.loadUi(ui_Path)

        fcSC_window.startBtn.clicked.connect(self.start_clicked)
        fcSC_window.quitBtn.clicked.connect(self.Deactivated)
        fcSC_window.resetBtn.clicked.connect(self.reset_clicked)
        fcSC_window.saveBtn.clicked.connect(self.save_clicked)
        fcSC_window.sumBtn.clicked.connect(self.sum_clicked)
        fcSC_window.totalRbtn.toggled.connect(self.btn_state)
        fcSC_window.incrRbtn.toggled.connect(self.btn_state)
        fcSC_window.averagedChk.toggled.connect(self.btn_state)
        fcSC_window.model1.toggled.connect(self.btn_state)
        fcSC_window.model2.toggled.connect(self.btn_state)

        fcSC_window.max_iter.textChanged.connect(self.max_iter_changed)
        fcSC_window.relax.textChanged.connect(self.relax_changed)
        fcSC_window.tableWidget.cellChanged.connect(self.cell_changed)

        fcSC_window.fy.setValidator(double_validator)
        fcSC_window.fc.setValidator(double_validator)
        fcSC_window.steps.setValidator(int_validator)
        fcSC_window.max_iter.setValidator(int_validator)
        fcSC_window.error.setValidator(double_validator)
        fcSC_window.relax.setValidator(double_validator)
        fcSC_window.target_LF_V.setValidator(double_validator)
        fcSC_window.target_LF_P.setValidator(double_validator)

        self.fc_default = "23.67"
        self.fy_default = "348.0"
        self.GZinput_default = "-10.0"
        self.rho_x_default = "0.005"
        self.rho_y_default = "0.005"
        self.rho_z_default = "0.005"
        self.steps_default = "10"
        self.max_iter_default = "20"
        self.error_default = "1.0e-2"
        self.relax_default = "1.2"
        self.scale_1_default = "2.0"
        self.scale_2_default = "1.2"
        self.scale_3_default = "1.2"
        self.target_LF_V_default = "2.0"
        self.target_LF_P_default = "1.0"
        self.disp_option_default = "incremental"
        self.eps_option_default = "eps_c"
        self.model_option_default = "1"
        self.averaged_Chk_default = "unaveraged"
        self.dx_default = "10.0"
        self.dy_default = "10.0"
        self.dz_default = "10.0"
        self.mat_values_default = [self.rho_x_default, self.rho_y_default, self.rho_z_default, self.dx_default,
                                   self.dy_default, self.dz_default]

        self.check_mat_objs()
        self.open_file()
        self.update_material()

        FCmw.addDockWidget(QtCore.Qt.RightDockWidgetArea, fcSC_window.dw)

    def check_mat_objs(self):
        #
        # mat_obj = {App::MaterialObjectPython : [rhox, rhoy, rhoz, dx, dy, dz]}
        #
        # objID = {row : App::MaterialObjectPython}
        #
        FreeCAD.ActiveDocument.recompute()
        self.objID = {}
        row = 0
        for obj in FreeCAD.ActiveDocument.Objects:
            if obj.TypeId == 'App::MaterialObjectPython' and obj.Label not in self.deleted_objs:
                self.objID[row] = obj
                row += 1

    def update_material(self):
        #
        # mat_obj = {App::MaterialObjectPython : [rhox, rhoy, rhoz, dx, dy, dz]}
        #
        # objID = {row : App::MaterialObjectPython}
        #
        from PySide2.QtWidgets import QTableWidgetItem
        from PySide2.QtGui import QFontMetrics
        max_width = 0
        tw = fcSC_window.tableWidget
        font_metrics = QFontMetrics(tw.font())
        tw.setRowCount(0)

        for row in self.objID:
            obj = self.objID[row]
            tw.setRowCount(row + 1)
            width = font_metrics.horizontalAdvance(obj.Label)
            if width > max_width: max_width = width
            if obj not in self.mat_obj:
                self.mat_obj[obj] = self.mat_values_default
            for field_index, field in enumerate(self.mat_obj[obj]):
                item = QTableWidgetItem(field)
                tw.setItem(row, field_index, item)

        tw.verticalHeader().setFixedWidth(max_width + 10)

        for row in self.objID:
            if obj.Label not in self.deleted_objs:
                tw.setVerticalHeaderItem(row, QTableWidgetItem(self.objID[row].Label))

    def cell_changed(self, r, c):
        self.mat_obj[self.objID[r]][c] = fcSC_window.tableWidget.item(r, c).text()

    def Deactivated(self):
        self.save_clicked()
        try:
            if fcSC_window.dw.isVisible():
                fcSC_window.dw.setVisible(False)
        except Exception:
            None

        FreeCAD.removeDocumentObserver(self.obs)

    def start_clicked(self):
        from PySide2 import QtWidgets

        self.save_clicked()

        try:
            rv = FCmw.findChild(QtWidgets.QTextEdit, "Report view")
            rv.clear()
        except Exception:
            None

        FreeCADGui.Selection.clearSelection()

        fcSC_macro = open(self.macro_file_path).read()
        exec(fcSC_macro)

    def quit_clicked(self):
        self.Deactivated()

    def reset_clicked(self):
        fcSC_window.max_iter.setText(self.max_iter_default)
        fcSC_window.error.setText(self.error_default)
        fcSC_window.relax.setText(self.relax_default)
        fcSC_window.incrRbtn.setChecked(True)
        fcSC_window.averagedChk.setChecked(False)

    def save_clicked(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '_sc.inp')

        with open(inp_file_path, "w") as f:
            f.write(fcSC_window.fc.text() + "\n")
            f.write(fcSC_window.fy.text() + "\n")
            f.write(fcSC_window.GZinput.text() + "\n")
            f.write(self.rho_x_default + "\n")
            f.write(self.rho_y_default + "\n")
            f.write(self.rho_z_default + "\n")
            f.write(fcSC_window.steps.text() + "\n")
            f.write(fcSC_window.max_iter.text() + "\n")
            f.write(fcSC_window.error.text() + "\n")
            f.write(fcSC_window.relax.text() + "\n")
            f.write(self.scale_1_default + "\n")
            f.write(self.scale_2_default + "\n")
            f.write(self.scale_3_default + "\n")
            f.write(self.disp_option + "\n")
            f.write(fcSC_window.target_LF_V.text() + "\n")
            f.write(self.eps_option + "\n")
            f.write(self.averaged_option + "\n")
            f.write(self.model_option + "\n")
            f.write(fcSC_window.target_LF_P.text() + "\n")
            f.write(str(len(self.mat_obj)) + "\n")
            for row in self.objID:
                for prop in self.mat_obj[self.objID[row]]:
                    f.write(prop + " ")
                f.write("\n")

    def sum_clicked(self):
        fcSC_sum = open(self.sum_file_path).read()
        exec(fcSC_sum)

        return

    def open_file(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '_sc.inp')

        # self.mat_obj = {key: self.mat_values_default for key in self.mat_obj}
        self.deleted_objs = []

        try:
            with open(inp_file_path, "r") as f:
                fcSC_window.fc.setText(str(f.readline().strip()))
                fcSC_window.fy.setText(str(f.readline().strip()))
                fcSC_window.GZinput.setText(str(f.readline().strip()))
                dummy = str(f.readline().strip())
                dummy = str(f.readline().strip())
                dummy = str(f.readline().strip())
                fcSC_window.steps.setText(str(f.readline().strip()))
                fcSC_window.max_iter.setText(str(f.readline().strip()))
                fcSC_window.error.setText(str(f.readline().strip()))
                fcSC_window.relax.setText(str(f.readline().strip()))
                dummy = str(f.readline().strip())
                dummy = str(f.readline().strip())
                dummy = str(f.readline().strip())
                if str(f.readline().strip()) == "total":
                    fcSC_window.totalRbtn.setChecked(True)
                else:
                    fcSC_window.incrRbtn.setChecked(True)
                LFinp_V = str(f.readline().strip())
                if LFinp_V == "":
                    fcSC_window.target_LF_V.setText(self.target_LF_V_default)
                else:
                    fcSC_window.target_LF_V.setText(LFinp_V)
                epsBtninp = f.readline().strip()
                avBtninp = str(f.readline().strip())
                mBtninp = str(f.readline().strip())
                if mBtninp == "1":
                    fcSC_window.model1.setChecked(True)
                if mBtninp == "2":
                    fcSC_window.model2.setChecked(True)
                LFinp_P = str(f.readline().strip())
                if LFinp_P == "":
                    fcSC_window.target_LF_P.setText(self.target_LF_P_default)
                else:
                    fcSC_window.target_LF_P.setText(LFinp_P)
                f.readline().strip()

                for row in self.objID:
                    values = str(f.readline().strip())
                    if values:
                        self.mat_obj[self.objID[row]] = values.split()
                    else:
                        self.mat_obj[self.objID[row]] = self.mat_values_default


        except FileNotFoundError:
            fcSC_window.fc.setText(self.fc_default)
            fcSC_window.fy.setText(self.fy_default)
            fcSC_window.GZinput.setText(self.GZinput_default)
            fcSC_window.steps.setText(self.steps_default)
            fcSC_window.max_iter.setText(self.max_iter_default)
            fcSC_window.error.setText(self.error_default)
            fcSC_window.relax.setText(self.relax_default)
            fcSC_window.incrRbtn.setChecked(True)
            fcSC_window.target_LF_V.setText(self.target_LF_V_default)
            fcSC_window.averagedChk.setChecked(False)
            fcSC_window.model1.setChecked(True)
            fcSC_window.target_LF_P.setText(self.target_LF_P_default)

    def max_iter_changed(self):
        if (fcSC_window.max_iter.text() != self.max_iter_default):
            fcSC_window.max_iter.setPalette(self.palette_warning)
        else:
            fcSC_window.max_iter.setPalette(self.palette_standard)

    def relax_changed(self):
        if (fcSC_window.relax.text() != self.relax_default):
            if fcSC_window.relax.text() == "":
                fcSC_window.relax.setText("0.0")
            if float(fcSC_window.relax.text()) > 1.5:
                fcSC_window.relax.setText("1.5")
            elif float(fcSC_window.relax.text()) < 1.0:
                fcSC_window.relax.setText("1.0")
            fcSC_window.relax.setPalette(self.palette_warning)
        else:
            fcSC_window.relax.setPalette(self.palette_standard)

    def btn_state(self):
        if fcSC_window.totalRbtn.isChecked():
            self.disp_option = "total"
        if fcSC_window.incrRbtn.isChecked():
            self.disp_option = "incremental"
        if fcSC_window.model1.isChecked():
            self.model_option = "1"
        if fcSC_window.model2.isChecked():
            self.model_option = "2"
        if fcSC_window.averagedChk.isChecked():
            self.averaged_option = "averaged"
        else:
            self.averaged_option = "unaveraged"


FreeCADGui.addWorkbench(fcSCWorkbench)
