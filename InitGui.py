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

        class DocObserver(object):  # document Observer
            def __init__(self, workbench_instance):
                self.workbench_instance = workbench_instance

            def slotActivateDocument(self, doc):
                print("slotActivateDocument")
                if FreeCAD.activeDocument().Label[0:7] != "Unnamed":
                    # print(self.workbench_instance.file_name)
                    self.workbench_instance.save_clicked()
                    self.workbench_instance.file_name = FreeCAD.activeDocument().Label
                    # print(self.workbench_instance.file_name)
                    self.workbench_instance.open_file()
                    fcSC_window.eps_c.setText("0.000")
                    fcSC_window.eps_s.setText("0.000")

            def slotFinishSaveDocument(self, doc, prop):
                print("slotFinishSaveDocument")
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
        fcSC_window.eps_cRbtn.toggled.connect(self.btn_state)
        fcSC_window.eps_sRbtn.toggled.connect(self.btn_state)
        fcSC_window.averagedChk.toggled.connect(self.btn_state)

        fcSC_window.max_iter.textChanged.connect(self.max_iter_changed)
        fcSC_window.relax.textChanged.connect(self.relax_changed)
        fcSC_window.scale_1.textChanged.connect(self.scale_1_changed)
        fcSC_window.scale_2.textChanged.connect(self.scale_2_changed)
        fcSC_window.scale_3.textChanged.connect(self.scale_3_changed)
        fcSC_window.rho_x.textChanged.connect(self.rho_x_changed)
        fcSC_window.rho_y.textChanged.connect(self.rho_y_changed)
        fcSC_window.rho_z.textChanged.connect(self.rho_z_changed)

        fcSC_window.fy.setValidator(double_validator)
        fcSC_window.fc.setValidator(double_validator)
        fcSC_window.rho_x.setValidator(double_validator)
        fcSC_window.rho_y.setValidator(double_validator)
        fcSC_window.rho_z.setValidator(double_validator)
        fcSC_window.steps.setValidator(int_validator)
        fcSC_window.max_iter.setValidator(int_validator)
        fcSC_window.error.setValidator(double_validator)
        fcSC_window.relax.setValidator(double_validator)
        fcSC_window.scale_1.setValidator(double_validator)
        fcSC_window.scale_2.setValidator(double_validator)
        fcSC_window.scale_3.setValidator(double_validator)
        fcSC_window.target_LF.setValidator(double_validator)

        self.fc_default = "23.67"
        self.fy_default = "348.0"
        self.GZinput_default = "-10.0"
        self.rho_x_default = "0.005"
        self.rho_y_default = "0.005"
        self.rho_z_default = "0.005"
        self.steps_default = "10"
        self.max_iter_default = "20"
        self.error_default = "1.0e-3"
        self.relax_default = "1.2"
        self.scale_1_default = "2.0"
        self.scale_2_default = "1.2"
        self.scale_3_default = "1.2"
        self.target_LF_default = "2.0"
        self.disp_option_default = "incremental"
        self.eps_option_default = "eps_c"
        self.averaged_Chk_default = "unaveraged"

        self.open_file()

        FCmw.addDockWidget(QtCore.Qt.RightDockWidgetArea, fcSC_window.dw)

    def Deactivated(self):

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

        print(self.macro_file_path)

        fcSC_macro = open(self.macro_file_path).read()
        exec(fcSC_macro)

    def quit_clicked(self):
        self.Deactivated()

    def reset_clicked(self):
        fcSC_window.max_iter.setText(self.max_iter_default)
        fcSC_window.error.setText(self.error_default)
        fcSC_window.relax.setText(self.relax_default)
        fcSC_window.scale_1.setText(self.scale_1_default)
        fcSC_window.scale_2.setText(self.scale_2_default)
        fcSC_window.scale_3.setText(self.scale_3_default)
        fcSC_window.relax.setPalette(self.palette_standard)
        fcSC_window.scale_1.setPalette(self.palette_standard)
        fcSC_window.scale_2.setPalette(self.palette_standard)
        fcSC_window.scale_3.setPalette(self.palette_standard)
        fcSC_window.incrRbtn.setChecked(True)
        fcSC_window.averagedChk.setChecked(False)

    def save_clicked(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '_sc.inp')
        with open(inp_file_path, "w") as f:
            f.write(fcSC_window.fc.text() + "\n")
            f.write(fcSC_window.fy.text() + "\n")
            f.write(fcSC_window.GZinput.text() + "\n")
            f.write(fcSC_window.rho_x.text() + "\n")
            f.write(fcSC_window.rho_y.text() + "\n")
            f.write(fcSC_window.rho_z.text() + "\n")
            f.write(fcSC_window.steps.text() + "\n")
            f.write(fcSC_window.max_iter.text() + "\n")
            f.write(fcSC_window.error.text() + "\n")
            f.write(fcSC_window.relax.text() + "\n")
            f.write(fcSC_window.scale_1.text() + "\n")
            f.write(fcSC_window.scale_2.text() + "\n")
            f.write(fcSC_window.scale_3.text() + "\n")
            f.write(self.disp_option + "\n")
            f.write(fcSC_window.target_LF.text() + "\n")
            f.write(self.eps_option + "\n")
            f.write(self.averaged_option + "\n")

    def sum_clicked(self):
        fcSC_sum = open(self.sum_file_path).read()
        exec(fcSC_sum)

        return

    def open_file(self):
        inp_file_path = os.path.join(self.dir_name, "control files", self.file_name + '_sc.inp')
        try:
            with open(inp_file_path, "r") as f:
                fcSC_window.fc.setText(str(f.readline().strip()))
                fcSC_window.fy.setText(str(f.readline().strip()))
                fcSC_window.GZinput.setText(str(f.readline().strip()))
                fcSC_window.rho_x.setText(str(f.readline().strip()))
                fcSC_window.rho_y.setText(str(f.readline().strip()))
                fcSC_window.rho_z.setText(str(f.readline().strip()))
                fcSC_window.steps.setText(str(f.readline().strip()))
                fcSC_window.max_iter.setText(str(f.readline().strip()))
                fcSC_window.error.setText(str(f.readline().strip()))
                fcSC_window.relax.setText(str(f.readline().strip()))
                fcSC_window.scale_1.setText(str(f.readline().strip()))
                fcSC_window.scale_2.setText(str(f.readline().strip()))
                fcSC_window.scale_3.setText(str(f.readline().strip()))
                if str(f.readline().strip()) == "total":
                    fcSC_window.totalRbtn.setChecked(True)
                else:
                    fcSC_window.incrRbtn.setChecked(True)
                LFinp = str(f.readline().strip())
                if LFinp == "":
                    fcSC_window.target_LF.setText(self.target_LF_default)
                else:
                    fcSC_window.target_LF.setText(LFinp)
                epsBtninp = f.readline().strip()
                if epsBtninp == "eps_c":
                    fcSC_window.eps_cRbtn.setChecked(True)
                else:
                    fcSC_window.eps_sRbtn.setChecked(True)
                avBtninp = str(f.readline().strip())
                if avBtninp == "averaged":
                    fcSC_window.averagedChk.setChecked(True)
                else:
                    fcSC_window.averagedChk.setChecked(False)


        except FileNotFoundError:
            fcSC_window.fc.setText(self.fc_default)
            fcSC_window.fy.setText(self.fy_default)
            fcSC_window.GZinput.setText(self.GZinput_default)
            fcSC_window.rho_x.setText(self.rho_x_default)
            fcSC_window.rho_y.setText(self.rho_y_default)
            fcSC_window.rho_z.setText(self.rho_z_default)
            fcSC_window.steps.setText(self.steps_default)
            fcSC_window.max_iter.setText(self.max_iter_default)
            fcSC_window.error.setText(self.error_default)
            fcSC_window.relax.setText(self.relax_default)
            fcSC_window.scale_1.setText(self.scale_1_default)
            fcSC_window.scale_2.setText(self.scale_2_default)
            fcSC_window.scale_3.setText(self.scale_3_default)
            fcSC_window.incrRbtn.setChecked(True)
            fcSC_window.target_LF.setText(self.target_LF_default)
            fcSC_window.eps_cRbtn.setChecked(True)
            fcSC_window.averagedChk.setChecked(False)

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

    def scale_1_changed(self):
        if (fcSC_window.scale_1.text() != self.scale_1_default):
            if fcSC_window.scale_1.text() == "":
                fcSC_window.scale_1.setText("0.0")
            if float(fcSC_window.scale_1.text()) > 3.0:
                fcSC_window.scale_1.setText("3.0")
            elif float(fcSC_window.scale_1.text()) < 1.0:
                fcSC_window.scale_1.setText("1.0")
            fcSC_window.scale_1.setPalette(self.palette_warning)
        else:
            fcSC_window.scale_1.setPalette(self.palette_standard)

    def scale_2_changed(self):
        if (fcSC_window.scale_2.text() != self.scale_2_default):
            if fcSC_window.scale_2.text() == "":
                fcSC_window.scale_2.setText("0.0")
            if float(fcSC_window.scale_2.text()) > 2.0:
                fcSC_window.scale_2.setText("2.0")
            elif float(fcSC_window.scale_2.text()) < 1.0:
                fcSC_window.scale_2.setText("1.0")
            fcSC_window.scale_2.setPalette(self.palette_warning)
        else:
            fcSC_window.scale_2.setPalette(self.palette_standard)

    def scale_3_changed(self):
        if fcSC_window.scale_3.text() == "":
            fcSC_window.scale_3.setText("0.0")
        if (fcSC_window.scale_3.text() != self.scale_3_default):
            if float(fcSC_window.scale_3.text()) > 2.0:
                fcSC_window.scale_3.setText("2.0")
            elif float(fcSC_window.scale_3.text()) < 1.0:
                fcSC_window.scale_3.setText("1.0")
            fcSC_window.scale_3.setPalette(self.palette_warning)
        else:
            fcSC_window.scale_3.setPalette(self.palette_standard)

    def rho_x_changed(self):
        try:
            rx = float(fcSC_window.rho_x.text())
        except ValueError:
            rx = 0.0
        if rx < 0.0:
            fcSC_window.rho_x.setText(self.rho_x_default)
            rx = float(self.rho_x_default)
        if rx < float(self.rho_x_default):
            fcSC_window.rho_x.setPalette(self.palette_warning)
        else:
            fcSC_window.rho_x.setPalette(self.palette_standard)

    def rho_y_changed(self):
        try:
            ry = float(fcSC_window.rho_y.text())
        except ValueError:
            ry = 0.0
        if ry < 0.0:
            fcSC_window.rho_y.setText(self.rho_y_default)
            ry = float(self.rho_y_default)
        if ry < float(self.rho_y_default):
            fcSC_window.rho_y.setPalette(self.palette_warning)
        else:
            fcSC_window.rho_y.setPalette(self.palette_standard)

    def rho_z_changed(self):
        try:
            rz = float(fcSC_window.rho_z.text())
        except ValueError:
            rz = 0.0
        if rz < 0.0:
            fcSC_window.rho_z.setText(self.rho_z_default)
            rz = float(self.rho_z_default)
        if rz < float(self.rho_z_default):
            fcSC_window.rho_z.setPalette(self.palette_warning)
        else:
            fcSC_window.rho_z.setPalette(self.palette_standard)


    def btn_state(self):
        if fcSC_window.totalRbtn.isChecked():
            self.disp_option = "total"
        if fcSC_window.incrRbtn.isChecked():
            self.disp_option = "incremental"
        if fcSC_window.eps_cRbtn.isChecked():
            self.eps_option = "eps_c"
        if fcSC_window.eps_sRbtn.isChecked():
            self.eps_option = "eps_s"
        if fcSC_window.averagedChk.isChecked():
            self.averaged_option = "averaged"
        else:
            self.averaged_option = "unaveraged"


FreeCADGui.addWorkbench(fcSCWorkbench)
