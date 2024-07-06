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

import os
import time
import math
import dummySC
import FemGui
import FreeCAD
import FreeCADGui
import ObjectsFem
import numpy as np
import Part as Part
import pyvista as pv
import FreeCAD as App
import FreeCADGui as Gui
from FreeCAD import Units
import scipy.sparse as scsp
from numba import jit, types
from numba.typed import Dict
from pyvista import CellType
import matplotlib.pyplot as plt
from femtools import membertools
from femmesh import meshsetsgetter
from femmesh import meshtools as mt
from femresult import resulttools as rt
from feminout import importToolsFem as itf
from matplotlib.widgets import Button, TextBox
from matplotlib.ticker import FormatStrFormatter
from femtaskpanels import task_result_mechanical as trm

global mdir
mdir = os.path.dirname(dummySC.file_path())

print("fcSC.py")

settings = {}
try:
    with open(os.path.join(mdir, 'fcVM.ini'), "r") as f:
        key = str(f.readline().strip()).split(" #")[0]
        settings[key] = int(f.readline().strip())
except FileNotFoundError:
    print("File fcVM.ini not found")

# print("settings: ", settings)

if settings["solver"] == 1:
    from sksparse.cholmod import cholesky
elif settings["solver"] == 2:
    from cholespy import CholeskySolverD, MatrixType
elif settings["solver"] == 3:
    from sksparse_minimal import SparseCholesky

global name
name = App.ActiveDocument.Label
file_path = os.path.join(mdir, "control files", name + '.inp')
macro_path = os.path.join(mdir, 'source code')
np.set_printoptions(precision=5, linewidth=300)


def prn_upd(*args):
    for obj in args:
        print(str(obj), end='')
    print()
    Gui.updateGui()


def setUpAnalysis():
    doc = App.ActiveDocument
    try:
        mesh = doc.getObject("FEMMeshGmsh").FemMesh
    except:
        mesh = None
    return_code = 0
    if mesh is None:
        # prn_upd("No Gmsh object. Please create one first")
        return_code = 1
        # raise SystemExit()
    elif mesh.Nodes == {}:
        # prn_upd("No mesh. Please generate mesh first")
        return_code = 2
        # raise SystemExit()

    analysis = doc.getObject("Analysis")
    if analysis is None:
        # prn_upd("No Analysis object. Please create one first")
        return_code = 3
        # raise SystemExit()

    # purge result objects
    rt.purge_results(analysis)
    doc.recompute()

    return doc, mesh, analysis, return_code


def setUpInput(doc, mesh, analysis):
    analysis = doc.getObject("Analysis")
    solver = doc.getObject("SolverCcxTools")
    docmesh = doc.getObject("FEMMeshGmsh")
    member = membertools.AnalysisMember(analysis)

    if solver == None:
        FemGui.setActiveAnalysis(App.activeDocument().Analysis)
        FemGui.getActiveAnalysis().addObject(ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument))
        solver = doc.getObject("SolverCcxTools")

    # determine elements connected to a node using FC API
    fet = mt.get_femelement_table(mesh)
    # fet is dictionary: { elementid : [ nodeid, nodeid, ... , nodeid ] }
    net = mt.get_femnodes_ele_table(mesh.Nodes, fet)
    # net is dictionary: {nodeID : [[eleID, binary node position], [], ...], nodeID : [[], [], ...], ...}
    # node0 has binary node position 2^0 = 1, node1 = 2^1 = 2, ..., node10 = 2^10 = 1024

    # create connectivity array elNodes for mapping local node number -> global node number
    elNodes = np.array([mesh.getElementNodes(el) for el in mesh.Volumes])  # elNodes[elementIndex] = [node1,...,Node10]

    # create nodal coordinate array nocoord for node number -> (x,y,z)
    ncv = list(mesh.Nodes.values())
    nocoord = np.asarray([[v.x, v.y, v.z] for v in ncv])  # nocoord[nodeIndex] = [x-coord, y-coord, z-coord]

    # get access to element sets: meshdatagetter.mat_geo_sets
    meshdatagetter = meshsetsgetter.MeshSetsGetter(
        analysis,
        solver,
        docmesh,
        member)
    meshdatagetter.get_mesh_sets()

    if len(member.mats_linear) == 1:
        element_sets = [mesh.Volumes]
    else:
        element_sets = [es["FEMElements"] for es in member.mats_linear]

    matCon = {}  # BooleanFragment Primitive the material object refers to
    ppEl = {}  # BooleanFragment Primitive element El belongs to
    materialbyElement = []  # see further

    prn_upd("Number of material objects: ", len(member.mats_linear))

    for indm, matobject in enumerate(member.mats_linear):
        E = float(App.Units.Quantity(matobject['Object'].Material['YoungsModulus']).getValueAs('MPa'))
        Nu = float(matobject['Object'].Material['PoissonRatio'])
        Density = float(App.Units.Quantity(matobject['Object'].Material['Density']).getValueAs('kg/mm^3'))
        prn_upd("Material Object: ", matobject['Object'].Name, "   E= ", E, "   Nu= ", Nu, "   Density= ", Density)
        for el in element_sets[indm]:  # element_sets[indm]: all elements with material indm
            if matCon: ppEl[el] = matCon[indm]  # ppEl[el]: primitive el belongs to
            materialbyElement.append([E, Nu, Density])  # materialbyElement[elementIndex] = [E, Nu, Density]

    materialbyElement = np.asarray(materialbyElement)

    noce = np.zeros((len(nocoord)), dtype=np.int16)
    for i in net:
        noce[i - 1] = len(net[i])

    # create boundary condition array dispfaces
    u0 = Units.Quantity(0.0, 1)
    dispfaces = []
    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintFixed') or obj.isDerivedFrom('Fem::ConstraintDisplacement'):

            bcnodes = []

            if obj.isDerivedFrom('Fem::ConstraintFixed'):
                bctype = [False, False, False]
                bcvalue = [u0, u0, u0]
            else:
                bctype = [obj.xFree, obj.yFree, obj.zFree]
                bcvalue = [obj.xDisplacement, obj.yDisplacement, obj.zDisplacement]

            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        bc = mesh.getNodesByVertex(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Edge:
                        bc = mesh.getNodesByEdge(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Face:
                        bc = mesh.getNodesByFace(
                            ref)  # all nodes on a primitive face with a displacement boundary condition
                        for bcn in bc:
                            bcnodes.append(bcn)
                    else:
                        prn_upd("No Boundaries Found")
            bcnodes = list(dict.fromkeys(bcnodes))  # remove duplicates in bcnodes
            if bcnodes:
                dispfaces.append([bcnodes, bctype, bcvalue])

    fix = Dict.empty(key_type=types.int64, value_type=types.float64)
    nn = 3 * len(nocoord)
    fixdof = np.ones(nn, dtype=int)
    movdof = np.zeros(nn, dtype=int)

    for face in dispfaces:
        if not face[1][0]:
            for node in face[0]:
                dof = 3 * (node - 1)
                try:
                    val = face[2][0].Value
                except:
                    val = face[2][0]
                fix[dof] = val
                fixdof[dof] = 0
        if not face[1][1]:
            for node in face[0]:
                dof = 3 * (node - 1) + 1
                try:
                    val = face[2][1].Value
                except:
                    val = face[2][1]
                fix[dof] = val
                fixdof[dof] = 0
        if not face[1][2]:
            for node in face[0]:
                dof = 3 * (node - 1) + 2
                try:
                    val = face[2][2].Value
                except:
                    val = face[2][2]
                fix[dof] = val
                fixdof[dof] = 0

    for dof in fix:
        if fix[dof] != 0.0:
            movdof[dof] = 1

    lf = [[0, 0, 0, 0, 0, 0]]  # load face nodes - signature for numba
    pr = [0.0]  # load face pressure - signature for numba
    lf_vertex = [[0]]
    pr_vertex = [[0.0, 0.0, 0.0]]
    lf_edge = [[0, 0, 0]]
    pr_edge = [[0.0, 0.0, 0.0]]
    lf_face = [[0, 0, 0, 0, 0, 0]]
    pr_face = [[0.0, 0.0, 0.0]]

    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintPressure'):
            if obj.Reversed:
                sign = 1
            else:
                sign = -1
            for part, faces in obj.References:  # obj.References: references to loaded primitive faces
                for face in faces:
                    ref = part.Shape.getElement(face)
                    if type(ref) == Part.Face:
                        for faceID in mesh.getFacesByFace(ref):  # face ID: ID of a 6-node face element
                            face_nodes = list(mesh.getElementNodes(faceID))  # 6-node element node numbers
                            lf.append(face_nodes)
                            pr.append(sign * obj.Pressure)
                    else:
                        prn_upd("No Faces with Pressure Loads")

        if obj.isDerivedFrom('Fem::ConstraintForce'):
            F = obj.Force
            d = obj.DirectionVector
            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        dp = [F * d.x, F * d.y, F * d.z]
                        lf_vertex.append(list(mesh.getNodesByVertex(ref)))
                        pr_vertex.append(dp)
                    elif type(ref) == Part.Edge:
                        L = ref.Length
                        dl = [F * d.x / L, F * d.y / L, F * d.z / L]
                        for edgeID in mesh.getEdgesByEdge(ref):
                            lf_edge.append(list(mesh.getElementNodes(edgeID)))
                            pr_edge.append(dl)
                    elif type(ref) == Part.Face:
                        A = ref.Area
                        dp = [F * d.x / A, F * d.y / A, F * d.z / A]
                        for faceID in mesh.getFacesByFace(ref):
                            lf_face.append(list(mesh.getElementNodes(faceID)))
                            pr_face.append(dp)
                    else:
                        prn_upd("No Boundaries Found")

    loadfaces = np.array(lf)
    pressure = np.array(pr)
    loadvertices = np.array(lf_vertex)
    vertexloads = np.array(pr_vertex)
    loadedges = np.array(lf_edge)
    edgeloads = np.array(pr_edge)
    loadfaces_uni = np.array(lf_face)
    faceloads = np.array(pr_face)

    # re-order element nodes
    for el in elNodes:
        el[1], el[2] = el[2], el[1]
        el[4], el[6] = el[6], el[4]
        el[8], el[9] = el[9], el[8]

    a_c = np.zeros(4 * len(elNodes), dtype=np.uint8)
    a_s = np.ones(12 * len(elNodes), dtype=np.float64)
    ev = np.zeros(24 * len(elNodes), dtype=np.float64)

    return (elNodes, nocoord, fix, fixdof, movdof, materialbyElement, noce,
            loadfaces, pressure,
            loadvertices, vertexloads,
            loadedges, edgeloads,
            loadfaces_uni, faceloads, a_c, a_s, ev)


# shape functions for a 4-node tetrahedron - only used for stress interpolation
@jit(nopython=True, cache=True)
def shape4tet(xi, et, ze, xl):
    shp = np.zeros((4), dtype=np.float64)

    # shape functions
    shp[0] = 1.0 - xi - et - ze
    shp[1] = xi
    shp[2] = et
    shp[3] = ze

    return shp


@jit(nopython=True, cache=True)
def shp10tet(xi, et, ze):
    shp = np.zeros((10), dtype=np.float64)

    # shape functions - source: Calculix, G Dhondt
    a = 1.0 - xi - et - ze
    shp[0] = (2.0 * a - 1.0) * a
    shp[1] = xi * (2.0 * xi - 1.0)
    shp[2] = et * (2.0 * et - 1.0)
    shp[3] = ze * (2.0 * ze - 1.0)
    shp[4] = 4.0 * xi * a
    shp[5] = 4.0 * xi * et
    shp[6] = 4.0 * et * a
    shp[7] = 4.0 * ze * a
    shp[8] = 4.0 * xi * ze
    shp[9] = 4.0 * et * ze
    return shp


@jit(nopython=True, cache=True, fastmath=True)
def dshp10tet(xi, et, ze, xl):
    dshp = np.zeros((3, 10), dtype=np.float64)
    dshpg = np.zeros((3, 10), dtype=np.float64)
    xs = np.zeros((3, 3), dtype=np.float64)
    xsi = np.zeros((3, 3), dtype=np.float64)
    bmat = np.zeros((6, 30), dtype=np.float64)

    # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
    dshp[0][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    dshp[0][1] = 4.0 * xi - 1.0
    # dshp[0][2] = 0.0
    # dshp[0][3] = 0.0
    dshp[0][4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
    dshp[0][5] = 4.0 * et
    dshp[0][6] = -4.0 * et
    dshp[0][7] = -4.0 * ze
    dshp[0][8] = 4.0 * ze
    # dshp[0][9] = 0.0

    # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
    dshp[1][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    # dshp[1][1] = 0.0
    dshp[1][2] = 4.0 * et - 1.0
    # dshp[1][3] = 0.0
    dshp[1][4] = -4.0 * xi
    dshp[1][5] = 4.0 * xi
    dshp[1][6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
    dshp[1][7] = -4.0 * ze
    # dshp[1][8] = 0.0
    dshp[1][9] = 4.0 * ze

    # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
    dshp[2][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    # dshp[2][1] = 0.0
    # dshp[2][2] = 0.0
    dshp[2][3] = 4.0 * ze - 1.0
    dshp[2][4] = -4.0 * xi
    # dshp[2][5] = 0.0
    dshp[2][6] = -4.0 * et
    dshp[2][7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
    dshp[2][8] = 4.0 * xi
    dshp[2][9] = 4.0 * et

    # xs = np.dot(xl, dshp.T) # local derivative of the global coordinates

    for i in range(3):
        for j in range(3):
            xs[i][j] = 0.0
            for k in range(10):
                xs[i][j] += xl[k][i] * dshp[j][k]

    # xsj = np.linalg.det(xs) # Jacobian

    xsj = (xs[0][0] * xs[1][1] * xs[2][2] -
           xs[0][0] * xs[1][2] * xs[2][1] +
           xs[0][2] * xs[1][0] * xs[2][1] -
           xs[0][2] * xs[1][1] * xs[2][0] +
           xs[0][1] * xs[1][2] * xs[2][0] -
           xs[0][1] * xs[1][0] * xs[2][2])

    # xsi = np.linalg.inv(xs) # global derivative of the local coordinates

    xsi[0][0] = (xs[1][1] * xs[2][2] - xs[2][1] * xs[1][2]) / xsj
    xsi[0][1] = (xs[0][2] * xs[2][1] - xs[0][1] * xs[2][2]) / xsj
    xsi[0][2] = (xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1]) / xsj
    xsi[1][0] = (xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2]) / xsj
    xsi[1][1] = (xs[0][0] * xs[2][2] - xs[0][2] * xs[2][0]) / xsj
    xsi[1][2] = (xs[1][0] * xs[0][2] - xs[0][0] * xs[1][2]) / xsj
    xsi[2][0] = (xs[1][0] * xs[2][1] - xs[2][0] * xs[1][1]) / xsj
    xsi[2][1] = (xs[2][0] * xs[0][1] - xs[0][0] * xs[2][1]) / xsj
    xsi[2][2] = (xs[0][0] * xs[1][1] - xs[1][0] * xs[0][1]) / xsj

    # dshp = np.dot(xsi.T, dshp) # global derivatives of the shape functions

    for i in range(3):
        for j in range(10):
            for k in range(3):
                dshpg[i][j] += xsi[k][i] * dshp[k][j]

    # computation of the strain interpolation matrix bmat
    for i in range(10):
        i3 = 3 * i
        d00 = dshpg[0][i]
        d10 = dshpg[1][i]
        d20 = dshpg[2][i]
        bmat[0][i3] = d00
        bmat[1][i3 + 1] = d10
        bmat[2][i3 + 2] = d20
        bmat[3][i3] = d10
        bmat[3][i3 + 1] = d00
        bmat[4][i3] = d20
        bmat[4][i3 + 2] = d00
        bmat[5][i3 + 1] = d20
        bmat[5][i3 + 2] = d10

    return xsj, bmat


# shape functions and their derivatives for a 6-node triangular interface element
@jit(nopython=True, cache=True)
def shape6tri(xi, et, xl):
    shp = np.zeros((6), dtype=np.float64)
    dshp = np.zeros((2, 6), dtype=np.float64)
    bmat = np.zeros((3, 36), dtype=np.float64)

    # shape functions
    shp[0] = (1.0 - xi - et) * (1.0 - 2.0 * xi - 2.0 * et)
    shp[1] = xi * (2.0 * xi - 1.0)
    shp[2] = et * (2.0 * et - 1.0)
    shp[3] = 4.0 * xi * (1.0 - xi - et)
    shp[4] = 4.0 * xi * et
    shp[5] = 4.0 * et * (1 - xi - et)

    # local derivatives of the shape functions: xi-derivative
    dshp[0][0] = -3.0 + 4.0 * et + 4.0 * xi
    dshp[0][1] = -1.0 + 4.0 * xi
    dshp[0][2] = 0.0
    dshp[0][3] = -4.0 * (-1.0 + et + 2.0 * xi)
    dshp[0][4] = 4.0 * et
    dshp[0][5] = -4.0 * et

    # local derivatives of the shape functions: eta-derivative
    dshp[1][0] = -3.0 + 4.0 * et + 4.0 * xi
    dshp[1][1] = 0.0
    dshp[1][2] = -1.0 + 4.0 * et
    dshp[1][3] = -4.0 * xi
    dshp[1][4] = 4.0 * xi
    dshp[1][5] = -4.0 * (-1.0 + 2.0 * et + xi)

    xs = np.dot(dshp, xl.T)  # xs = [ [[dx/dxi],[dy/dxi],[dz/dxi]] , [[dx/det],[dy/det],[dz/det]] ]

    xp = np.cross(xs[0], xs[1])  # vector normal to surface

    xsj = np.linalg.norm(xp)  # Jacobian

    # xsj = np.sqrt(xp[0]*xp[0]+xp[1]*xp[1]+xp[2]*xp[2])

    xx = xs[0] / np.linalg.norm(xs[0])  # unit vector in xi direction

    # xx = np.sqrt(xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2])

    xp /= xsj  # unit vector normal to surface
    xt = np.cross(xp, xx)  # unit vector tangential to surface and normal to xx

    # computation of the "strain" interpolation matrix bmat
    for i in range(6):
        ia = 3 * i
        ib = ia + 18
        ni = shp[i]
        bmat[0][ia] = ni
        bmat[1][ia + 1] = ni
        bmat[2][ia + 2] = ni
        bmat[0][ib] = -ni
        bmat[1][ib + 1] = -ni
        bmat[2][ib + 2] = -ni

    return xsj, shp, bmat, xx, xt, xp


@jit(nopython=True, cache=True)
def shape2lin(xi, xle):
    shp = np.zeros((3), dtype=np.float64)
    dshp = np.zeros((3), dtype=np.float64)

    # shape functions
    shp[0] = - 0.5 * (1.0 - xi) * xi
    shp[1] = 0.5 * (1.0 + xi) * xi
    shp[2] = (1.0 + xi) * (1.0 - xi)

    # local derivatives of the shape functions: xi-derivative
    dshp[0] = xi - 0.5
    dshp[1] = xi + 0.5
    dshp[2] = - 2.0 * xi

    dx_dxi = xle[0][0] * dshp[0] + xle[0][1] * dshp[1] + xle[0][2] * dshp[2]
    dy_dxi = xle[1][0] * dshp[0] + xle[1][1] * dshp[1] + xle[1][2] * dshp[2]
    dz_dxi = xle[2][0] * dshp[0] + xle[2][1] * dshp[1] + xle[2][2] * dshp[2]

    xsj = np.sqrt(dx_dxi ** 2 + dy_dxi ** 2 + dz_dxi ** 2)

    return xsj, shp


# linear-elastic material stiffness matrix

@jit(nopython=True, cache=True)
def hooke(element, materialbyElement, active_c, active_s, ev, rhox, rhoy, rhoz, Es):
    margin = 1.0
    dmat = np.zeros((6, 6), dtype=np.float64)
    # e = materialbyElement[element][0] + (rhox + rhoy + rhoz) * Es / 3
    # e = materialbyElement[element][0] + max(rhox, rhoy, rhoz) * Es / 3
    e = materialbyElement[element][0]
    dmat[0][0] = dmat[1][1] = dmat[2][2] = e
    dmat[3][3] = dmat[4][4] = dmat[5][5] = e / 2.0

    u = ev[0:3]
    v = ev[3:6]

    df1ds = np.array([u[0] ** 2, u[1] ** 2, u[2] ** 2, 2 * u[0] * u[1], 2 * u[2] * u[0], 2 * u[1] * u[2]])
    df2ds = np.array([v[0] ** 2, v[1] ** 2, v[2] ** 2, 2 * v[0] * v[1], 2 * v[2] * v[0], 2 * v[1] * v[2]])
    # df1ds = np.array([u[0] ** 2, u[1] ** 2, u[2] ** 2, u[0] * u[1], u[2] * u[0], u[1] * u[2]])
    # df2ds = np.array([v[0] ** 2, v[1] ** 2, v[2] ** 2, v[0] * v[1], v[2] * v[0], v[1] * v[2]])

    # df1ds = np.array([u[0] ** 2, u[1] ** 2, u[2] ** 2, u[0] * u[1], u[1] * u[2], u[2] * u[0]])
    # print("np.dot(dmat, df1ds)", np.dot(dmat, df1ds))

    # print("u: ", u)
    # print("v: ", v)
    # print("active_c: ", active_c)

    if active_c == 1:  # one yield surface
        # print("df1ds: ", df1ds)
        dfdf = np.dot(df1ds.reshape(1, 6).T, df1ds.reshape(1, 6))
        # print("dfdf: ")
        # for row in dfdf:
        #     print(row)
        d = np.dot(df1ds, np.dot(dmat, df1ds))
        M = np.eye(6) - np.dot(dfdf, dmat) / (margin * d)
        # dmat = dmat - e * dfdf
        # print("M: ")
        # for row in M:
        #     print(row)
        # print("e, d: ", e, d)
        dmat = np.dot(dmat, M)
        # print("dmat: ")
        # for row in dmat:
        #     print(row)


    elif active_c == 2:  # two yield surfaces
        df1df1 = np.dot(df1ds.reshape(1, 6).T, df1ds.reshape(1, 6))
        df2df2 = np.dot(df2ds.reshape(1, 6).T, df2ds.reshape(1, 6))
        df1df2 = np.dot(df1ds.reshape(1, 6).T, df2ds.reshape(1, 6))
        df2df1 = np.dot(df2ds.reshape(1, 6).T, df1ds.reshape(1, 6))

        a11 = np.dot(df1ds, np.dot(dmat, df1ds))
        a12 = np.dot(df1ds, np.dot(dmat, df2ds))
        a21 = np.dot(df2ds, np.dot(dmat, df1ds))
        a22 = np.dot(df2ds, np.dot(dmat, df2ds))

        d = a11 * a22 - a12 * a21

        dfdf = a22 * df1df1 + a11 * df2df2 - a12 * df1df2 - a21 * df2df1

        M = np.eye(6) - np.dot(dfdf, dmat) / (margin * d)
        dmat = np.dot(dmat, M)
    elif active_c == 3:  # apex
        dmat = (1.0 - 1.0 / margin) * dmat

    # print(active_c)

    # print("active_c: ", active_c)

    # dmat[0][0] += active_s[0] * rhox * Es
    # dmat[1][1] += active_s[1] * rhoy * Es
    # dmat[2][2] += active_s[2] * rhoz * Es

    # dmat[0][0] +=  rhox * Es
    # dmat[1][1] +=  rhoy * Es
    # dmat[2][2] +=  rhoz * Es

    return dmat


# Gaussian integration points and weights
@jit(nopython=True, cache=True)
def gaussPoints():
    # Gaussian integration points and weights for 10-noded tetrahedron
    gp10 = np.array([[0.138196601125011, 0.138196601125011, 0.138196601125011,
                      0.041666666666667],
                     [0.585410196624968, 0.138196601125011, 0.138196601125011,
                      0.041666666666667],
                     [0.138196601125011, 0.585410196624968, 0.138196601125011,
                      0.041666666666667],
                     [0.138196601125011, 0.138196601125011, 0.585410196624968,
                      0.041666666666667]])
    # Gaussian integration points and weights for 6-noded triangle
    gp6 = np.array([[0.445948490915965, 0.445948490915965,
                     0.111690794839005],
                    [0.10810301816807, 0.445948490915965,
                     0.111690794839005],
                    [0.445948490915965, 0.10810301816807,
                     0.111690794839005],
                    [0.091576213509771, 0.091576213509771,
                     0.054975871827661],
                    [0.816847572980458, 0.091576213509771,
                     0.054975871827661],
                    [0.091576213509771, 0.816847572980458,
                     0.054975871827661]])
    # Gaussian integration points and weights for 3-noded line
    gp2 = np.array([[-0.5773502691896257, 1.0], [0.5773502691896257, 1.0]])

    return gp10, gp6, gp2


# calculate the global stiffness matrix and load vector
# @jit(
#     "types.Tuple((float64[:],int64[:], int64[:], float64[:], float64[:]))(int64[:,:], float64[:,:], float64[:,:], DictType(int64,float64), int64[:,:], float64, float64, float64, float64[:])",
#     nopython=True, cache=True)
@jit(nopython=True, cache=True)
def calcGSM(elNodes, nocoord, materialbyElement, fix, grav_z, rhox, rhoy, rhoz, loadfaces, pressure,
            loadvertices, vertexloads, loadedges, edgeloads, loadfaces_uni, faceloads, a_c, a_s, ev, model):
    gp10, gp6, gp2 = gaussPoints()
    ne = len(elNodes)  # number of volume elements
    # nn = len(nocoord[:, 0])  # number of degrees of freedom
    nn = len(nocoord)  # number of degrees of freedom

    # number of entries in the lower diagonal of the stiffness matrix: (dof**2 + dof) / 2
    ns = int((30 * 30 + 30) / 2 * ne)

    print(f"ne = {ne}")
    print(f"ns = {ns}\n")

    xle = np.zeros((3, 3), dtype=types.float64)  # coordinates of load line nodes
    xlf = np.zeros((3, 6), dtype=types.float64)  # coordinates of load face nodes
    xlv = np.zeros((10, 3), dtype=types.float64)  # coordinates of volume element nodes
    glv = np.zeros((3 * nn), dtype=types.float64)  # global load vector
    row = np.zeros(ns, dtype=types.int64)  # row indices of COO matrix
    col = np.zeros(ns, dtype=types.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=types.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=types.float64)  # modification to the stiffness matrix for displacement BCs
    dof = np.zeros(30, dtype=types.int64)
    bmatV = np.zeros((6, 30), dtype=np.float64)
    x = np.zeros((4 * ne, 3), dtype=np.float64)

    #   calculate element load vectors for pressure and add to global vector

    for face in range(len(pressure) - 1):  # first pressure value is a dummy signature for numba
        if len(pressure) == 1:
            break

        nda = loadfaces[face + 1]  # node numbers of loaded face
        for i in range(3):
            for j in range(6):
                nd = nda[j]
                xlf[i][j] = nocoord[nd - 1][i]  # coordinates of loaded face nodes

        # integrate element load vector
        for index in range(len(gp6)):
            xi = gp6[index][0]
            et = gp6[index][1]
            xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
            nl = 0
            for i in range(len(loadfaces[face] - 1)):
                nd = loadfaces[face + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                for k in range(3):
                    load = shp[nl] * pressure[face + 1] * xp[k] * abs(xsj) * gp6[index][2]
                    glv[iglob3 + k] += load
                nl += 1

    for vertex in range(len(loadvertices) - 1):  # first vertex is a dummy signature for numba
        if len(loadvertices) == 1:
            break

        nd = loadvertices[vertex + 1][0]
        iglob = nd - 1
        iglob3 = 3 * iglob
        glv[iglob3:iglob3 + 3] += vertexloads[vertex + 1]

    for face in range(len(loadfaces_uni) - 1):  # first face is a dummy signature for numba
        if len(loadfaces_uni) == 1:
            break

        nda = loadfaces_uni[face + 1]  # node numbers of loaded face
        for i in range(3):
            for j in range(6):
                nd = nda[j]
                xlf[i][j] = nocoord[nd - 1][i]  # coordinates of loaded face nodes

        # integrate element load vector
        for index in range(len(gp6)):
            xi = gp6[index][0]
            et = gp6[index][1]
            xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
            nl = 0
            for i in range(len(loadfaces_uni[face] - 1)):
                nd = loadfaces_uni[face + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * faceloads[face + 1] * abs(xsj) * gp6[index][2]
                glv[iglob3:iglob3 + 3] += load
                nl += 1

    for edge in range(len(loadedges) - 1):  # first edge is a dummy signature for numba
        if len(loadedges) == 1:
            break

        nda = loadedges[edge + 1]  # node numbers of loaded edge
        for i in range(3):
            for j in range(3):
                nd = nda[j]
                xle[i][j] = nocoord[nd - 1][i]  # coordinates of loaded edge nodes
        # integrate element load vector
        for index in range(len(gp2)):
            xi = gp2[index][0]
            xsj, shp = shape2lin(xi, xle)
            nl = 0
            for i in range(len(loadedges[edge] - 1)):
                nd = loadedges[edge + 1][i]
                iglob = nd - 1
                iglob3 = 3 * iglob
                load = shp[nl] * edgeloads[edge + 1] * abs(xsj) * gp2[index][1]
                glv[iglob3:iglob3 + 3] += load
                nl += 1

    # for each volume element calculate the element stiffness matrix
    # and gravity load vector and add to global matrix and vector

    pos = 0

    V = 0.0

    density = materialbyElement[0][2]

    for el, nodes in enumerate(elNodes):
        esm = np.zeros((30, 30), dtype=np.float64)
        gamma = np.zeros((30), dtype=np.float64)

        # set up nodal values for this element
        for i in range(3):
            for j in range(10):
                nd = nodes[j]
                xlv[j][i] = nocoord[nd - 1][i]

        # integrate element matrix
        for i, ip in enumerate(gp10):
            # dmat = np.zeros((6, 6), dtype=np.float64)
            #
            # dmat = hooke(0, materialbyElement, 0,
            #              np.array(3 * [1.0]), np.array(6 * [1.0]), rhox, rhoy, rhoz, 210000.0)

            dmatc = np.zeros((6, 6), dtype=np.float64)
            Ec = materialbyElement[el][0]
            dmatc[0][0] = dmatc[1][1] = dmatc[2][2] = Ec
            dmatc[3][3] = dmatc[4][4] = dmatc[5][5] = Ec / 2.0

            dmats = np.zeros((6, 6), dtype=np.float64)
            Es = 210000.0
            dmats[0, 0] = rhox * Es
            dmats[1, 1] = rhoy * Es
            dmats[2, 2] = rhoz * Es

            if model == 1:
                dmat = dmatc + dmats
            elif model == 2:
                dmat = dmatc + dmats

            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shp10tet(xi, et, ze)
            xsj, bmatV = dshp10tet(xi, et, ze, xlv)
            esm += np.dot(bmatV.T, np.dot(dmat, bmatV)) * ip[3] * abs(xsj)
            gamma[2::3] += grav_z * density * shp * ip[3] * abs(xsj)
            V += xsj * ip[3]  # Element volume - not used
            x[4 * el + i] = np.dot(xlv.T, shp)  #

        for i in range(10):
            nd = nodes[i] - 1
            glv[3 * nd] += gamma[3 * i]
            glv[3 * nd + 1] += gamma[3 * i + 1]
            glv[3 * nd + 2] += gamma[3 * i + 2]
            for j in range(3):
                dof[3 * i + j] = 3 * nd + j

        for i in range(30):
            dofi = dof[i]
            if dofi in fix:
                row[pos] = dofi
                col[pos] = dofi
                stm[pos] = 1.0
                modf[dofi] += fix[dofi]
                pos += 1
                for j in range(i):
                    dofj = dof[j]
                    if dofj not in fix:
                        modf[dofj] -= esm[i][j] * fix[dofi]
            else:
                for j in range(i + 1):
                    dofj = dof[j]
                    if dofj in fix:
                        modf[dofi] -= esm[i][j] * fix[dofj]
                    else:
                        if dofi > dofj:
                            row[pos] = dofi
                            col[pos] = dofj
                        else:
                            row[pos] = dofj
                            col[pos] = dofi
                        stm[pos] = esm[i][j]
                        pos += 1

    row_r = row[:pos]
    col_r = col[:pos]
    stm_r = stm[:pos]

    loadsumx = 0.0
    loadsumy = 0.0
    loadsumz = 0.0
    for node in range(nn):
        dof = 3 * node
        loadsumx += glv[dof]
        loadsumy += glv[dof + 1]
        loadsumz += glv[dof + 2]

    print("mesh volume", V)
    print("loadsumx", loadsumx)
    print("loadsumy", loadsumy)
    print("loadsumz", loadsumz, "\n")

    # for rw in dmat:
    #     print(rw)
    # print()

    # print("ns, pos: ", ns, pos)

    return stm_r, row_r, col_r, glv, modf, V, loadsumx, loadsumy, loadsumz, ne, nn, x


@jit(nopython=True, cache=True)
def calcTSM(elNodes, nocoord, materialbyElement, fix, rhox, rhoy, rhoz, a_c, a_s, ev):
    print("-------------------------calcTSM-------------------------")
    gp10, gp6, gp2 = gaussPoints()
    ne = len(elNodes)  # number of volume elements
    # nn = len(nocoord[:, 0])  # number of degrees of freedom
    nn = len(nocoord)  # number of nodes

    # number of entries in the lower diagonal of the stiffness matrix: (dof**2 + dof) / 2
    ns = int((30 * 30 + 30) / 2 * ne)

    # print(f"ne = {ne}")
    # print(f"ns = {ns}\n")

    print(a_s[0:3])

    xlv = np.zeros((10, 3), dtype=types.float64)  # coordinates of volume element nodes
    row = np.zeros(ns, dtype=types.int64)  # row indices of COO matrix
    col = np.zeros(ns, dtype=types.int64)  # column indices of COO matrix
    stm = np.zeros(ns, dtype=types.float64)  # stiffness values of COO matrix
    modf = np.zeros((3 * nn), dtype=types.float64)  # modification to the stiffness matrix for displacement BCs
    dof = np.zeros(30, dtype=types.int64)

    # for each volume element calculate the element stiffness matrix
    # and add to global matrix and vector

    pos = 0

    for el, nodes in enumerate(elNodes):
        esm = np.zeros((30, 30), dtype=np.float64)

        # set up nodal values for this element
        for i in range(3):
            for j in range(10):
                nd = nodes[j]
                xlv[j][i] = nocoord[nd - 1][i]

        # integrate element matrix
        for i, ip in enumerate(gp10):
            # dmat = np.zeros((6, 6), dtype=np.float64)
            # print(rhox, rhoy, rhoz)

            dmatc = hooke(0, materialbyElement, a_c[4 * el + i],
                          a_s[12 * el + 3 * i:12 * el + 3 * i + 3], ev[24 * el + 6 * i:24 * el + 6 * i + 6], rhox, rhoy,
                          rhoz,
                          210000.0)

            dmats = np.zeros((6, 6), dtype=np.float64)
            Es = 210000.0
            dmats[0, 0] = rhox * Es
            dmats[1, 1] = rhoy * Es
            dmats[2, 2] = rhoz * Es

            # dmat = hooke(0, materialbyElement, 0,
            #              np.array(3 * [1.0]), np.array(6 * [1.0]), rhox, rhoy, rhoz, 210000.0)

            # print(a_c[4 * el + i])
            # print(dmat[0][0])

            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            xsj, bmatV = dshp10tet(xi, et, ze, xlv)
            esm += np.dot(bmatV.T, np.dot(dmatc + dmats, bmatV)) * ip[3] * abs(xsj)

        for i in range(10):
            nd = nodes[i] - 1
            for j in range(3):
                dof[3 * i + j] = 3 * nd + j

        for i in range(30):
            dofi = dof[i]
            if dofi in fix:
                row[pos] = dofi
                col[pos] = dofi
                stm[pos] = 1.0
                modf[dofi] += fix[dofi]
                pos += 1
                for j in range(i):
                    dofj = dof[j]
                    if dofj not in fix:
                        modf[dofj] -= esm[i][j] * fix[dofi]
            else:
                for j in range(i + 1):
                    dofj = dof[j]
                    if dofj in fix:
                        modf[dofi] -= esm[i][j] * fix[dofj]
                    else:
                        if dofi > dofj:
                            row[pos] = dofi
                            col[pos] = dofj
                        else:
                            row[pos] = dofj
                            col[pos] = dofi
                        stm[pos] = esm[i][j]
                        pos += 1

    row_r = row[:pos]
    col_r = col[:pos]
    stm_r = stm[:pos]

    print("dmatc")
    for rw in dmatc:
        print(rw)
    print()

    # print("esm")
    # for rw in esm:
    #     print(rw)
    #     print()
    # print()
    #
    #
    # print("ns, pos: ", ns, pos)

    return stm_r, row_r, col_r, modf


# calculate load-deflection curve
def calcDisp(elNodes, nocoord, fixdof, movdof, modf, materialbyElement, stm, row, col,
             glv, nstep, iterat_max, error_max, relax, scale_re, scale_up, scale_dn, fcd, sig_yield_inp, disp_output,
             fcSC_window, target_LF, x, noce, rhox, rhoy, rhoz, a_c, a_s, ev, fix, model):
    ndof = len(glv)  # number of degrees of freedom
    nelem = len(elNodes)  # number of elements

    pr_u = fcSC_window.progressBar.setValue
    st_u = fcSC_window.Step.setText
    LF_u = fcSC_window.Load_Factor.setText
    EC_u = fcSC_window.eps_c.setText
    ES_u = fcSC_window.eps_s.setText

    t0 = time.perf_counter()
    gsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
    t1 = time.perf_counter()
    prn_upd("construct sparse global stiffness matrix: {:.2e} s".format((t1 - t0)))

    # TODO: improve error estimate for pure displacement control
    qnorm = np.linalg.norm(glv)
    # if qnorm < 1.0: qnorm = 1.0

    # Cholesky decomposition of the global stiffness matrix and elastic solution using Cholmod
    t0 = time.perf_counter()
    if settings["solver"] == 1:
        factor = cholesky(gsm)
    elif settings["solver"] == 2:
        solver = CholeskySolverD(gsm.shape[0], gsm.indptr, gsm.indices, gsm.data, MatrixType.CSC)
    elif settings["solver"] == 3:
        sparse_cholesky = SparseCholesky(gsm)

    t1 = time.perf_counter()
    f = fixdof * glv + modf
    if settings["solver"] == 1:
        ue = factor(f)  # elastic solution
    elif settings["solver"] == 2:
        ue = np.empty(gsm.shape[0])
        solver.solve(f, ue)
    elif settings["solver"] == 3:
        ue = sparse_cholesky.solve_A(f)

    # print("ue_max (after GSM): ", np.max(np.absolute(ue)))

    # gsmd = gsm.todense().A
    # gsmd = np.tril(gsmd) + np.tril(gsmd, -1).T
    # Q = np.dot(gsmd, ue)
    # print("ERROR Elastic Stiffness Matrix: ", np.linalg.norm(Q - f) / np.linalg.norm(f))

    t2 = time.perf_counter()
    prn_upd("sparse Cholesky decomposition: {:.2e} s, elastic solution: {:.2e} s".format((t1 - t0), (t2 - t1)))

    # initiate analysis
    dl0 = min(target_LF, 1.0) / nstep
    if float(nstep) == 1.0 : dl0 = target_LF # nstep == 1 execute an elastic analysis
    dl = dl0
    du = dl * ue

    sig_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    sig_old = np.zeros(24 * nelem, dtype=np.float64)
    sigc_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    sigc_old = np.zeros(24 * nelem, dtype=np.float64)
    sigs_new = np.zeros(24 * nelem, dtype=np.float64)  # stress in Tet10
    sigs_old = np.zeros(24 * nelem, dtype=np.float64)
    epscc_old = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic compressive strain concrete
    epscc_new = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic compressive strain concrete
    epsct_old = np.zeros(4 * nelem, dtype=np.float64)
    epsct_new = np.zeros(4 * nelem, dtype=np.float64)
    epss_new = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic compressive strain concrete
    epss_old = np.zeros(4 * nelem, dtype=np.float64)
    sig_yield = np.full(4 * nelem, sig_yield_inp, dtype=np.float64)  # yield stress in Tet10
    sig_test = np.zeros(24 * nelem, dtype=np.float64)
    peeq = np.zeros(4 * nelem, dtype=np.float64)  # equivalent plastic strain in Tet10
    triax = np.zeros(4 * nelem, dtype=np.float64)  # triaxiality in Tet10
    pressure = np.zeros(4 * nelem, dtype=np.float64)  # pressure in Tet10
    sigmises = np.zeros(4 * nelem, dtype=np.float64)  # von Mises stress in Tet10
    ecr = np.zeros(4 * nelem, dtype=np.float64)  # critical plastic strain in Tet10
    csr = np.zeros(4 * nelem, dtype=np.float64)  # critical strain ratio in Tet10
    disp_new = np.zeros(ndof, dtype=np.float64)  # displacement results
    disp_old = np.zeros(ndof, dtype=np.float64)  # displacement results
    lbd = np.zeros(1, dtype=np.float64)  # load level
    rfl = np.zeros(1, dtype=np.float64)  # reaction force level (for displacement control)
    qin = np.zeros(3 * len(nocoord), dtype=np.float64)

    gp10, gp6, gp2 = gaussPoints()

    # determine elastic reaction force on moving boundary
    if max(movdof) == 1:
        qelastic = np.zeros(3 * len(nocoord), dtype=np.float64)
        update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, ue, sig_old, sig_new,
                           sigc_old, sigc_new, sigs_old, sigs_new, epscc_old, epsct_old, epscc_new, epsct_new, epss_old, epss_new, sig_test,
                           qelastic, rhox, rhoy, rhoz, a_c, a_s, ev, 0, 0, model)

        qelastic *= movdof
        qelnorm = np.linalg.norm(qelastic)
        qnorm = qelnorm
        sig_new = np.zeros(24 * nelem, dtype=np.float64)  # reset sig_new to zero
        sigc_new = np.zeros(24 * nelem, dtype=np.float64)  # reset sigc_new to zero
        sigs_new = np.zeros(24 * nelem, dtype=np.float64)  # reset sigs_new to zero
        epscc_new = np.zeros(24 * nelem, dtype=np.float64)  # reset epscc_new to zero
        epsct_new = np.zeros(24 * nelem, dtype=np.float64)  # reset epsct_new to zero
        epss_new = np.zeros(24 * nelem, dtype=np.float64)  # reset epss_new to zero

    step = -1
    cnt = True
    fail = False
    dof_max = np.argmax(np.abs(ue))

    un = [0.]
    epsccplot = [0.]
    epsctplot = [0.]
    epssplot = [0.]
    crip = [0]
    pplot = [0.]
    svmplot = [0.]
    pdfcdplot = [0.]
    peeqplot = [0.]
    ecrplot = [0.]
    lout = [0.]

    if float(nstep) == 1.0:
        # perform an elastic (one-step) analysis
        step = 0
        out_disp = 1
        lbd = np.append(lbd, 1.0)
        rfl = np.append(rfl, 1.0)
        disp_new = dl * ue
        un.append(np.max(np.abs(disp_new)))
        qin = np.zeros(3 * len(nocoord), dtype=np.float64)
        sig_yield *= 1.0e6
        fcd *= 1.0e6
        update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, disp_new, sig_old, sig_new,
                           sigc_old, sigc_new, sigs_old, sigs_new, epscc_old, epsct_old, epscc_new, epsct_new, epss_old, epss_new, sig_test,
                           qin, rhox, rhoy, rhoz, a_c, a_s, ev, 0, 0, model)

        cnt = False

    factor_time_tot = 0.0
    iterat_tot = 0

    # fex_sum = 0.0
    # for node in range(len(nocoord)):
    #     dof = 3 * node
    #     fex_sum += f[dof]
    # print("(fixdof * glv + modf)_sum (after GSM): ", fex_sum)
    # print("ue_max (after GSM): ", np.max(ue))
    # print("elastic stiffness GSM: ", fex_sum / np.max(ue))

    error = 0.0

    while cnt:
        cnt = False
        iRiks = True
        pstep = 0
        pr_u(0)
        while pstep < nstep:
            step += 1
            pstep += 1
            st_u(str(pstep))
            restart = 0

            prn_upd("\n***************************** Step: {} *****************************".format(step))
            if iRiks:
                a = du.copy()  # a: Riks control vector
                # print("np.linalg.norm(a)", np.linalg.norm(a))
                a /= np.linalg.norm(a)
                sig_old = sig_new.copy()
                sigc_old = sigc_new.copy()
                sigs_old = sigs_new.copy()
                epscc_old = epscc_new.copy()
                epsct_old = epsct_new.copy()
                epss_old = epss_new.copy()
                lbd = np.append(lbd, lbd[step] + dl)  # lbd: load level
            else:
                lbd[step + 1] = lbd[step] + dl

            # update stresses and loads

            qin_old = qin.copy()

            qin = np.zeros(3 * len(nocoord), dtype=np.float64)

            iterat = 0

            # print("du_0_max: ", np.max(du))

            update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, du, sig_old, sig_new,
                               sigc_old, sigc_new, sigs_old, sigs_new, epscc_old, epsct_old, epscc_new, epsct_new, epss_old, epss_new, sig_test,
                               qin, rhox, rhoy, rhoz, a_c, a_s, ev,
                               step, iterat, model)

            # print("sig_new_max: ", np.max(sig_new))

            fex = fixdof * lbd[step + 1] * glv
            fin = fixdof * qin
            r = fex - fin

            rnorm = np.linalg.norm(r)
            # qnorm = np.linalg.norm(fin)

            # out-of-balance error

            error_old = error

            error = rnorm / qnorm

            prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))

            # tangent stiffness matrix

            # if error > error_max and model == 2:
            #     stm, row, col, modf = calcTSM(elNodes, nocoord, materialbyElement, fix, rhox, rhoy, rhoz, a_c,
            #                                   a_s, ev)
            #     t0 = time.perf_counter()
            #     tsm = scsp.csc_matrix((stm, (row, col)), shape=(ndof, ndof))  # construct sparse global stiffness matrix
            #     t1 = time.perf_counter()
            #     prn_upd("construct sparse tangent stiffness matrix: {:.2e} s".format((t1 - t0)))
            #     factor = cholesky(tsm)
            #     f = fixdof * glv + modf
            #     ue = factor(f)  # tangent elastic solution

            # print("ue: ",ue)
            # end tangent stiffness matrix

            # print("|f| after TSM: ", np.linalg.norm(f))
            # print("|ue|: ", np.linalg.norm(ue))

            fex_sum = 0.0
            fey_sum = 0.0
            fez_sum = 0.0
            for node in range(len(nocoord)):
                dof = 3 * node
                fex_sum += f[dof]
                fey_sum += f[dof + 1]
                fez_sum += f[dof + 2]
            # print("fex_sum (after TSM): ", fex_sum)
            # print("fey_sum (after TSM): ", fey_sum)
            # print("fez_sum (after TSM): ", fez_sum)
            # print("ue_max (after TSM): ", np.max(np.absolute(ue)))
            # print("stiffness after TSM: ", fex_sum / np.max(ue))

            while error > error_max:

                # print(">>>> ---------------- displacement correction -------------------")

                iterat += 1
                iterat_tot += 1

                # relax = 0.9

                f = relax * r
                t0 = time.perf_counter()
                if settings["solver"] == 1:
                    due = factor(f)
                elif settings["solver"] == 2:
                    due = np.empty(gsm.shape[0])
                    solver.solve(f, due)
                elif settings["solver"] == 3:
                    due = sparse_cholesky.solve_A(f)

                # print("|r|: ", np.linalg.norm(r))

                t1 = time.perf_counter()

                factor_time_tot += t1 - t0

                # print("norm(due): ", np.linalg.norm(due))

                # Riks control correction to load level increment

                # print("dl before riks: ", dl)
                # print("before riks lbd[step + 1]: ", lbd[step + 1])

                if iRiks:
                    teller = -np.dot(a, due)
                    noemer = np.dot(a, ue)
                    # print(">>>> lbd[step + 1]: ", lbd[step + 1])
                    # print(">>>> teller: ", teller, "noemer: ", noemer)
                    # print(">>>> teller/noemer: ", teller / noemer)
                    # dl = -np.dot(a, due) / np.dot(a, ue)
                    dl = teller / noemer
                    # if error < 0.99 * error_old: dl = 0.0
                    lbd[step + 1] += dl
                    # print(">>>> lbd[step + 1]: ", lbd[step + 1])
                else:
                    dl = 0.0

                # print("dl after riks: ", dl)
                # print("lbd[step + 1]: ", lbd[step + 1])
                # Riks control correction to displacement increment
                du += due + dl * ue

                # r_sum = 0.0
                # for node in range(len(nocoord)):
                #     dof = 3 * node
                #     r_sum += f[dof]
                # print("r_sum: ", r_sum)
                # print("due_max: ", np.max(due))

                # print("norm(du): ", np.linalg.norm(du))

                # update stresses and loads

                qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                # print(">>>> ---------------- update stresses and loads -------------------")

                update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, du, sig_old, sig_new,
                                   sigc_old, sigc_new, sigs_old, sigs_new, epscc_old, epsct_old, epscc_new, epsct_new, epss_old, epss_new,
                                   sig_test, qin, rhox, rhoy, rhoz, a_c, a_s,
                                   ev, step, iterat, model)

                # print("norm(qin): ", np.linalg.norm(qin))

                # calculate out of balance error
                r = fixdof * (lbd[step + 1] * glv - qin)
                rnorm = np.linalg.norm(r)
                # qnorm = np.linalg.norm(fixdof * qin)
                error = rnorm / qnorm

                prn_upd("Iteration: {}, Error: {:.2e}".format(iterat, error))
                # print(">>>> fex_norm: ", np.linalg.norm(fex))
                # print(">>>> fex_min: ", np.min(fex))
                # print(">>>> fex_max: ", np.max(fex))
                # print(">>>> fin_norm: ", np.linalg.norm(fin))
                # print(">>>> lbd[step+1]: ", lbd[step + 1])
                # print(">>>> rnorm: ", rnorm)
                # print(">>>> qnorm: ", qnorm)

                if iterat > iterat_max:
                    # scale down
                    if restart == 4:
                        print("MAXIMUM RESTARTS REACHED")
                        fail = True
                        return disp_new, sig_new, epscc_new, epsct_new, sigmises, epss_new, lout, un, crip, peeqplot, pplot, svmplot, pdfcdplot, epsccplot, epsctplot, epssplot, fail, pressure / fcd

                    restart += 1
                    if step > 0:
                        dl = (lbd[step] - lbd[step - 1]) / scale_re / restart
                        du = (disp_new - disp_old) / scale_re / restart
                    else:
                        # for first step only
                        dl = dl0 / scale_re / restart
                        du = dl * ue / scale_re / restart
                    lbd[step + 1] = lbd[step] + dl

                    qin = np.zeros(3 * len(nocoord), dtype=np.float64)
                    update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, du, sig_old, sig_new,
                                       sigc_old, sigc_new, sigs_old, sigs_new, epscc_old, epsct_old, epscc_new, epsct_new, epss_old, epss_new,
                                       sig_test, qin, rhox, rhoy, rhoz, a_c,
                                       a_s, ev, step, iterat, model)

                    r = fixdof * (lbd[step + 1] * (glv + modf) - qin)
                    rnorm = np.linalg.norm(r)
                    # qnorm = np.linalg.norm(fixdof * qin)
                    error = rnorm / qnorm

                    iterat = 0

            # if abs(lbd[step + 1]) > abs(target_LF) and iRiks:
            if abs(target_LF - lbd[step]) < abs(lbd[step + 1] - lbd[step]) and iRiks:
                print("REACHED TARGET LOAD")
                iRiks = False
                dl = target_LF - lbd[step]
                du = dl / (lbd[step + 1] - lbd[step]) * du
                step -= 1
                pstep -= 1
            else:
                # update results at end of converged load step
                pr_u(int(100 * (pstep + 1) / nstep))
                LF_u(str(round(lbd[step + 1], 3)))
                disp_old = disp_new.copy()
                disp_new += du
                dl = lbd[step + 1] - lbd[step]
                if max(movdof) == 1:
                    rfl = np.append(rfl, np.sum(movdof * qin))

                if iterat > 10:
                    # scale down
                    dl /= scale_dn
                    du /= scale_dn
                if iterat < 5:
                    # scale up
                    dl *= scale_up
                    du *= scale_up

                # un.append(np.max(np.abs(disp_new)))
                un.append(disp_new[dof_max])
                # update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, us_s, peeq, csr,
                #                 triax,
                #                 pressure, sigmises, ecr)
                pressure = (sig_new[0::6] + sig_new[1::6] + sig_new[2::6]) / 3
                epsccplot.append(np.max(epscc_new))
                epsctplot.append(np.max(epsct_new))
                epssplot.append(np.max(epss_new))
                pdfcdplot.append(np.min(pressure / fcd))

                EC_u(str(round(max(epscc_new), 3)))
                ES_u(str(round(max(epss_new), 3)))

                if not iRiks: break

        if max(movdof) == 1:
            lout = rfl
        else:
            lout = lbd

        prn_upd("\n***************************** Progress Report Step: {} *****************************".format(step))

        print(
            '{0: >10}{1: >10}{2: >10}{3: >10}{4: >10}{5: >10}'.format(
                "load",
                "un",
                "p/fcd",
                "eps_cc",
                "eps_ct",
                "eps_s"))
        for i in range(len(lout)):
            print(
                '{0: >10.2e}{1: >10.2e}{2: >10.2e}{3: >10.2e}{4: >10.2e}{5: >10.2e}'.format(
                    lout[i],
                    un[i],
                    pdfcdplot[i],
                    epsccplot[i],
                    epsctplot[i],
                    epssplot[i]))

        us_c = 0.0035
        us_s = 0.025

        if fcSC_window.eps_cRbtn.isChecked():
            epsc_limit = np.argwhere(np.asarray(epsccplot) > us_c)
            epsc_non_zero = np.nonzero(epsccplot)
            if len(epsc_non_zero[0]) != 0:
                el_limit = epsc_non_zero[0][0] - 1
            else:
                el_limit = 0
            if len(epsc_limit) != 0:
                ul_limit = epsc_limit[0][0] - 1
            else:
                ul_limit = 0
        else:
            epss_limit = np.argwhere(np.asarray(epssplot) > us_s)
            epss_non_zero = np.nonzero(epssplot)
            if len(epss_non_zero[0]) != 0:
                el_limit = epss_non_zero[0][0] - 1
            else:
                el_limit = 0
            if len(epss_limit) != 0:
                ul_limit = epss_limit[0][0] - 1
            else:
                ul_limit = 0

        averaged = fcSC_window.averagedChk.isChecked()
        cnt, dl, du, target_LF = plot(fcSC_window, averaged, el_limit, ul_limit, un, lout, epsccplot, epsctplot, epssplot, dl, du,
                                      target_LF, nstep, ue, us_s, disp_new, disp_old, elNodes, nocoord, sig_new,
                                      epss_new, sigmises, epscc_new, epsct_new, noce, sig_yield_inp, pressure / fcd)

        # print("target_LF: ", target_LF)
        # print("dl: ", dl)

    prn_upd("total time evaluating K_inv * r: {}".format(factor_time_tot))
    prn_upd("total number of iterations: {}".format(iterat_tot))
    if iterat_tot != 0:
        prn_upd("average time to evaluate K_inv * r: {}".format(factor_time_tot / iterat_tot))

    out = step + 1
    u_out = un[out] - un[out - 1]
    l_out = lbd[out] - lbd[out - 1]
    prn_upd("Step: {0:2d} Load level increment: {1:.3f} Displacement increment: {2:.4e}".format(out, l_out,
                                                                                                u_out))

    if disp_output == "total":
        return disp_new, sig_new, epscc_new, epsct_new, sigmises, epss_new, lout, un, crip, peeqplot, pplot, svmplot, pdfcdplot, epsccplot, epsctplot, epssplot, fail, pressure / fcd
    else:
        return disp_new - disp_old, sig_new, epscc_new, epsct_new, sigmises, epss_new, lout, un, crip, peeqplot, pplot, svmplot, pdfcdplot, epsccplot, epsctplot, epssplot, fail, pressure / fcd


# plot the load-deflection curve
def plot(fcVM, averaged, el_limit, ul_limit, un, lbd, epsccplot, epsctplot, epssplot, dl, du, target_LF, nstep, ue, us_s,
         disp_new, disp_old, elNodes, nocoord, sig_new, epss, sigmises, epscc, epsct, noce, sig_yield, pdfcd):

    class Index(object):

        def __init__(self, averaged, disp_new,
                     disp_old, elNodes, nocoord, sig_new, epss,
                     sigmises, epscc, epsct, noce, sig_yield):
            self.averaged = averaged
            self.elNodes = elNodes
            self.nocoord = nocoord
            self.sig_new = sig_new
            self.epss = epss
            self.sigmises = sigmises
            self.epscc = epscc
            self.epsct = epsct
            self.noce = noce
            self.disp_old = disp_old
            self.disp_new = disp_new
            self.sig_yield = sig_yield
            self.pdfcd = pdfcd

        def stop(self, event):
            self.cnt = False
            plt.close()
            self.clicked = True

        def add(self, event):
            self.cnt = True
            self.clicked = True
            print("self.LF: ", self.LF)
            print("self.target_LF: ", self.target_LF)
            print("self.target_LF_out: ", self.target_LF_out)
            cr1 = (self.target_LF - self.LF) * (self.target_LF_out - self.LF) <= 0.0
            print(cr1)
            if (cr1):
                self.dl = np.sign(self.target_LF_out - self.LF) * 1.0 / self.nstep
                self.du = self.dl * self.ue
                print("self.dl: ", self.dl)
            plt.close()

        def rev(self, event):
            self.cnt = True
            self.dl = - self.dl
            self.du = - self.du
            plt.close()
            self.clicked = True

        def close_window(self, event):
            if self.cnt == False:
                self.stop('stop_event')

        def submit(self, LF):
            self.target_LF_out = float(LF)

        def PSV(self, event):

            class Toggle_Plane():
                def __init__(self, p):
                    pass

                def __call__(self, *args, **kwargs):
                    p.plane_widgets[0].SetEnabled(not p.plane_widgets[0].GetEnabled())

                def update(self):
                    pass

            class Scale_Stress():
                def __init__(self, p, scale):
                    self.scale = scale

                def __call__(self, value):
                    glyphs1 = grid1.glyph(orient="Major Principal Stress Vector", scale="Major Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs1, name='sv1', show_scalar_bar=False, lighting=False, cmap=['red'])

                    glyphs2 = grid1.glyph(orient="Intermediate Principal Stress Vector",
                                          scale="Intermediate Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs2, name='sv2', show_scalar_bar=False, lighting=False, cmap=['green'])

                    glyphs3 = grid1.glyph(orient="Minor Principal Stress Vector", scale="Minor Principal Stress",
                                          factor=10 ** value * self.scale,
                                          geom=geom)

                    p.add_mesh(glyphs3, name='sv3', show_scalar_bar=False, lighting=False, cmap=['blue'])

                def update(self):
                    pass

            def screen_shot(p):
                file = os.path.join(mdir, "output files", name + '_PSV.png')
                p.screenshot(file)

            # pv.global_theme.cmap = 'jet'
            pv.global_theme.cmap = 'coolwarm'
            # pv.global_theme.cmap = 'turbo'
            # pv.set_plot_theme('dark')
            pv.global_theme.full_screen = True
            pv.global_theme.title = 'VTK'

            # Controlling the text properties
            sargs = dict(
                title_font_size=16,
                label_font_size=14,
                shadow=False,
                color="white",
                bold=True,
                n_labels=5,
                italic=True,
                fmt="%.2e",
                font_family="arial",

            )

            tet10stress, tet10epss, tet10epscc, tet10epsct, tet10svm, tet10pdfcd = mapStresses(self.averaged, self.elNodes,
                                                                                  self.nocoord,
                                                                                  self.sig_new, self.epscc, self.epsct,
                                                                                  self.sigmises, self.epss, self.noce,
                                                                                  self.sig_yield, self.pdfcd)

            tet10s1, tet10s2, tet10s3, sv1, sv2, sv3 = calculate_principal_stress(tet10stress)

            x_range = max(nocoord[:, 0]) - min(nocoord[:, 0])
            y_range = max(nocoord[:, 1]) - min(nocoord[:, 1])
            z_range = max(nocoord[:, 2]) - min(nocoord[:, 2])

            geo_range = max(x_range, y_range, z_range)

            stress_range = max(max(map(abs, tet10s1)), max(map(abs, tet10s2)), max(map(abs, tet10s3)))

            if stress_range == 0.0:
                scale = 1.0
            else:
                scale = 0.2 * geo_range / stress_range

            disp_range = max(self.disp_new) - min(self.disp_new)

            if disp_range == 0.0:
                scale_disp = 1.0
            else:
                scale_disp = 0.2 * geo_range / disp_range

            padding = np.full(len(self.elNodes), 10, dtype=int)
            self.elements = np.vstack((padding, (self.elNodes - 1).T)).T

            celltypes = np.full(len(self.elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)
            geom = pv.Line()

            points = self.nocoord

            grid = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid.point_data['Displacement'] = np.reshape(disp_new, (len(nocoord), 3))

            grid1 = grid.warp_by_vector(vectors='Displacement', factor=scale_disp,
                                        inplace=False, progress_bar=False)
            grid1.point_data["Major Principal Stress"] = tet10s1.flatten(order="F")
            grid1.point_data["Intermediate Principal Stress"] = tet10s2.flatten(order="F")
            grid1.point_data["Minor Principal Stress"] = tet10s3.flatten(order="F")
            grid1.point_data['Major Principal Stress Vector'] = sv1
            grid1.point_data['Intermediate Principal Stress Vector'] = sv2
            grid1.point_data['Minor Principal Stress Vector'] = sv3

            glyphs1 = grid1.glyph(orient="Major Principal Stress Vector", scale="Major Principal Stress", factor=scale,
                                  geom=geom)
            glyphs2 = grid1.glyph(orient="Intermediate Principal Stress Vector", scale="Intermediate Principal Stress",
                                  factor=scale, geom=geom)
            glyphs3 = grid1.glyph(orient="Minor Principal Stress Vector", scale="Minor Principal Stress", factor=scale,
                                  geom=geom)

            p = pv.Plotter()

            p.set_background('grey', all_renderers=True)

            # p.set_background("royalblue", top="aliceblue")

            p.add_mesh_clip_plane(grid1, name="mesh", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  show_scalar_bar=False,
                                  scalar_bar_args=sargs, cmap=['grey'])
            p.add_mesh(glyphs1, name='sv1', show_scalar_bar=False, lighting=False, cmap=['red'])
            p.add_mesh(glyphs2, name='sv2', show_scalar_bar=False, lighting=False, cmap=['green'])
            p.add_mesh(glyphs3, name='sv3', show_scalar_bar=False, lighting=False, cmap=['blue'])

            p.add_key_event("s", lambda: screen_shot(p))

            ss = Scale_Stress(p, scale)
            p.add_slider_widget(
                callback=ss,
                rng=[-1.0, 1.0],
                value=0.0,
                title="log(Scale Stress)",
                pointa=(0.05, 0.075),
                pointb=(0.3, 0.075),
                style='modern',
                slider_width=0.02,
                tube_width=0.02
            )

            tp = Toggle_Plane(p)
            p.add_key_event("t", lambda: tp())

            p.show(cpos=[1.0, 1.0, 1.0])

        def VTK(self, event):

            class Toggle_Plane():
                def __init__(self, p):
                    self.normals = []
                    self.origins = []
                    self.names = ["mesh1", "mesh2", "mesh3", "mesh4"]

                def __call__(self, *args, **kwargs):
                    if p.plane_widgets:
                        self.normals = []
                        self.origins = []
                        for widget in p.plane_widgets:
                            self.normals.append(widget.GetNormal())
                            self.origins.append(widget.GetOrigin())
                        p.clear_plane_widgets()

                    else:
                        for name in self.names:
                            p.remove_actor(name)
                        p.subplot(0, 0)
                        p.add_mesh_clip_plane(grid1, name="mesh1", show_edges=False, normal=self.normals[0],
                                              origin=self.origins[0],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(0, 1)
                        p.add_mesh_clip_plane(grid2, name="mesh2", show_edges=False, normal=self.normals[1],
                                              origin=self.origins[1],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(1, 0)
                        p.add_mesh_clip_plane(grid3, name="mesh3", show_edges=False, normal=self.normals[2],
                                              origin=self.origins[2],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)
                        p.subplot(1, 1)
                        p.add_mesh_clip_plane(grid4, name="mesh4", show_edges=False, normal=self.normals[3],
                                              origin=self.origins[3],
                                              invert=True, lighting=True,
                                              scalar_bar_args=sargs)

                def update(self):
                    pass

            def screen_shot(p):
                file = os.path.join(mdir, "output files", name + '.png')
                p.screenshot(file)

            # pv.global_theme.cmap = 'jet'
            pv.global_theme.cmap = 'coolwarm'
            # pv.global_theme.cmap = 'turbo'
            # pv.set_plot_theme('dark')
            pv.global_theme.full_screen = True
            pv.global_theme.title = 'VTK'

            # Controlling the text properties
            sargs = dict(
                title_font_size=16,
                label_font_size=14,
                shadow=False,
                color="white",
                bold=True,
                n_labels=5,
                italic=True,
                fmt="%.2e",
                font_family="arial",

            )

            tet10stress, tet10epss, tet10epscc, tet10epsct, tet10svm, tet10pdfcd = mapStresses(self.averaged, self.elNodes,
                                                                                  self.nocoord,
                                                                                  self.sig_new, self.epscc, self.epsct,
                                                                                  self.sigmises, self.epss, self.noce,
                                                                                  self.sig_yield, self.pdfcd)

            x_range = max(nocoord[:, 0]) - min(nocoord[:, 0])
            y_range = max(nocoord[:, 1]) - min(nocoord[:, 1])
            z_range = max(nocoord[:, 2]) - min(nocoord[:, 2])

            geo_range = max(x_range, y_range, z_range)

            disp_range = max(self.disp_new) - min(self.disp_new)

            if disp_range == 0.0:
                scale = 1.0
            else:
                scale = 0.2 * geo_range / disp_range

            padding = np.full(len(self.elNodes), 10, dtype=int)
            self.elements = np.vstack((padding, (self.elNodes - 1).T)).T

            celltypes = np.full(len(self.elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)

            points = self.nocoord + scale * np.reshape(disp_new, (len(nocoord), 3))
            grid1 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid1.point_data["Concrete Plastic Compressive Strain\n"] = tet10epscc.flatten(order="F")
            grid2 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid2.point_data["Steel Plastic Tensile Strain\n"] = tet10epss.flatten(order="F")
            grid3 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid3.point_data["Concrete Plastic Tensile Strain\n"] = tet10epsct.flatten(order="F")
            grid4 = pv.UnstructuredGrid(self.elements, celltypes, points)
            grid4.point_data["p/fcd\n"] = tet10pdfcd.flatten(order="F")

            p = pv.Plotter(shape=(2, 2))

            p.link_views()

            p.set_background('grey', all_renderers=True)

            # left upper pane
            p.subplot(0, 0)
            p.add_mesh_clip_plane(grid1, name="mesh1", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # right upper pane
            p.subplot(0, 1)
            p.add_mesh_clip_plane(grid2, name="mesh2", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # left lower pane
            p.subplot(1, 0)
            p.add_mesh_clip_plane(grid3, name="mesh3", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            # right lower pane
            p.subplot(1, 1)
            p.add_mesh_clip_plane(grid4, name="mesh4", show_edges=False, normal=[1.0, 0., 0.], invert=True,
                                  scalar_bar_args=sargs)

            p.add_key_event("s", lambda: screen_shot(p))

            tp = Toggle_Plane(p)
            p.add_key_event("t", lambda: tp())

            p.show(cpos=[1.0, 1.0, 1.0])

    print(len(un), len(lbd))
    callback = Index(averaged, disp_new,
                     disp_old, elNodes, nocoord, sig_new, epss,
                     sigmises, epscc, epsct, noce, sig_yield)
    callback.cnt = False
    callback.clicked = False
    callback.dl = dl
    callback.du = du
    callback.target_LF = target_LF
    callback.target_LF_out = target_LF
    callback.LF = lbd[-1]
    callback.ue = ue
    callback.nstep = nstep
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    fig.canvas.manager.set_window_title('fcSC')
    plt.subplots_adjust(bottom=0.2)
    ax[0].plot(un, lbd, '-ok')
    ax[0].set(xlabel='displacement [mm]', ylabel='load factor [-] or load [N]', title='')
    if fcVM.eps_cRbtn.isChecked():
        ax[1].set(xlabel='concrete plastic compressive strain [-]', title='')
        ax[1].plot(epsccplot, lbd, '-ok')
    else:
        ax[1].set(xlabel='steel plastic tensile strain [-]', title='')
        ax[1].plot(epssplot, lbd, '-ok')
    ax[0].grid()
    ax[1].grid()
    b_w = 0.075
    b_h = 0.06
    b_s = 0.01
    b_y = 0.05
    # axstop = plt.axes([0.5 - b_w - b_s, b_y, b_w, b_h])
    # axadd = plt.axes([0.5 + b_s, b_y, b_w, b_h])
    axstop = plt.axes([0.4 - b_w / 2.0 - b_w - b_s, b_y, b_w, b_h])
    axadd = plt.axes([0.4 - b_w / 2.0, b_y, b_w, b_h])
    axbox = plt.axes([0.4 - b_w / 2.0 + b_w + b_s, b_y, b_w, b_h])
    axVTK = plt.axes([0.675, b_y, b_w, b_h])
    axPSV = plt.axes([0.675 + b_w + b_s, b_y, b_w, b_h])
    # axbox = plt.axes([0.53, b_y, b_w, b_h])
    bstop = Button(axstop, 'stop')
    bstop.on_clicked(callback.stop)
    badd = Button(axadd, 'add')
    badd.on_clicked(callback.add)
    bVTK = Button(axVTK, 'VTK')
    bPSV = Button(axPSV, 'PSV')
    bVTK.on_clicked(callback.VTK)
    bPSV.on_clicked(callback.PSV)
    # brev = Button(axrev, 'rev')
    # brev.on_clicked(callback.rev)
    text_box = TextBox(axbox, "", textalignment="center")
    text_box.set_val(target_LF)
    text_box.on_submit(callback.submit)
    fig.text(0.4 - b_w / 2.0 + 2.0 * (b_w + b_s) - 0.005, b_y + b_h / 2.0 - 0.01, 'Target Load Factor', fontsize=10)
    fig.canvas.mpl_connect('close_event', callback.close_window)

    if ul_limit != 0:
        if fcVM.eps_cRbtn.isChecked():
            epsc_ult = 0.0035
            fac = (epsc_ult - epsccplot[ul_limit]) / (epsccplot[ul_limit + 1] - epsccplot[ul_limit])
            lbd_limit = lbd[ul_limit] + fac * (lbd[ul_limit + 1] - lbd[ul_limit])
            un_limit = un[ul_limit] + fac * (un[ul_limit + 1] - un[ul_limit])
            ax[0].plot([0.0, un_limit], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un_limit, un_limit], [0.0, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un[el_limit], un[el_limit]], [0.0, lbd[el_limit]], color='b', linestyle="--")
            ax[0].plot([0.0, un[el_limit]], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")
            ax[1].plot([0.0, epsc_ult], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[1].plot([epsc_ult, epsc_ult], [0.0, lbd_limit], color='r', linestyle="--")
            ax[1].plot([0.0, epsc_ult], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")
        else:
            # print("ul_limit: ", ul_limit)
            fac = (us_s - epssplot[ul_limit]) / (epssplot[ul_limit + 1] - epssplot[ul_limit])
            lbd_limit = lbd[ul_limit] + fac * (lbd[ul_limit + 1] - lbd[ul_limit])
            un_limit = un[ul_limit] + fac * (un[ul_limit + 1] - un[ul_limit])
            ax[0].plot([0.0, un_limit], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un_limit, un_limit], [0.0, lbd_limit], color='r', linestyle="--")
            ax[0].plot([un[el_limit], un[el_limit]], [0.0, lbd[el_limit]], color='b', linestyle="--")
            ax[0].plot([0.0, un[el_limit]], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")
            ax[1].plot([0.0, us_s], [lbd_limit, lbd_limit], color='r', linestyle="--")
            ax[1].plot([us_s, us_s], [0.0, lbd_limit], color='r', linestyle="--")
            ax[1].plot([0.0, us_s], [lbd[el_limit], lbd[el_limit]], color='b', linestyle="--")

    plt.show()

    while True:
        plt.pause(0.01)
        if callback.clicked:
            break

    return callback.cnt, callback.dl, callback.du, callback.target_LF_out


# update PEEQ and CSR
@jit(nopython=True, cache=True, fastmath=True)
def update_PEEQ_CSR(nelem, materialbyElement, sig_test, sig_new, sig_yield, us_s, peeq, csr, triax, pressure,
                    sigmises, ecr):
    E = materialbyElement[0][0]  # Young's Modulus
    nu = materialbyElement[0][1]  # Poisson's Ratio
    nu = 0.0  # for simplified concrete analysis
    G = E / 2.0 / (1 + nu)  # shear modulus
    if us_s == 0.0:
        us_s = 1.0e12
    alpha = np.sqrt(np.e) * us_s  # stress triaxiality T = 1/3 for uniaxial test
    beta = 1.5

    for el in range(nelem):
        for ip in range(4):
            ipos1 = 4 * el + ip
            ipos2 = 24 * el + 6 * ip
            st0, st1, st2, st3, st4, st5 = sig_test[ipos2:ipos2 + 6]
            sn0, sn1, sn2, sn3, sn4, sn5 = sig_new[ipos2:ipos2 + 6]
            p_t = (st0 + st1 + st2) / 3.0
            p_n = (sn0 + sn1 + sn2) / 3.0
            st0 -= p_t
            st1 -= p_t
            st2 -= p_t
            sn0 -= p_n
            sn1 -= p_n
            sn2 -= p_n
            sig_mises_test = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                                     3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))
            sig_mises_new = np.sqrt(1.5 * (sn0 ** 2 + sn1 ** 2 + sn2 ** 2) +
                                    3.0 * (sn3 ** 2 + sn4 ** 2 + sn5 ** 2))

            DL = 0.0
            if (sig_mises_test > sig_yield[ipos1]):
                DL = (sig_mises_test - sig_yield[ipos1]) / (3.0 * G)
                peeq[ipos1] += DL

            # if sig_mises_new > 0.0:
            #     T = p_n / sig_mises_new
            # else:
            #     T = 1.0e12

            T = p_n / sig_yield[ipos1]  # experimental - verify against theory and tests

            pressure[ipos1] = p_n
            sigmises[ipos1] = sig_mises_new
            triax[ipos1] = T

            critical_strain = alpha * np.exp(-beta * T)

            if critical_strain < 1.0e-6:
                critical_strain = 1.0e-6

            ecr[ipos1] = critical_strain

            csr[ipos1] += DL / critical_strain


# update stresses and loads

@jit(nopython=True, cache=True, fastmath=True)
def update_stress_load(gp10, elNodes, nocoord, materialbyElement, fcd, sig_yield, du, sig, sig_update,
                       sigc, sigc_update, sigs, sigs_update, epscc, epsct, epscc_update, epsct_update, epss, epss_update, sig_test_global,
                       qin, rhox, rhoy, rhoz, a_c, a_s, ev, step, iterat, model):
    # print("-------------------------update_stress_load-------------------------")
    u10 = np.empty(30, dtype=np.float64)  # displacements for the 10 tetrahedral nodes
    sig_test = np.empty(6, dtype=np.float64)
    sigc_test = np.empty(6, dtype=np.float64)
    sigs_test = np.empty(6, dtype=np.float64)
    bm0 = np.empty(10, dtype=np.float64)
    bm1 = np.empty(10, dtype=np.float64)
    bm2 = np.empty(10, dtype=np.float64)
    bm3 = np.empty(10, dtype=np.float64)
    bm4 = np.empty(10, dtype=np.float64)
    bm5 = np.empty(10, dtype=np.float64)
    bm6 = np.empty(10, dtype=np.float64)
    bm7 = np.empty(10, dtype=np.float64)
    bm8 = np.empty(10, dtype=np.float64)
    xlv0 = np.empty(10, dtype=np.float64)
    xlv1 = np.empty(10, dtype=np.float64)
    xlv2 = np.empty(10, dtype=np.float64)
    dshp0 = np.empty(10, dtype=np.float64)
    dshp1 = np.empty(10, dtype=np.float64)
    dshp2 = np.empty(10, dtype=np.float64)
    dshpg0 = np.empty(10, dtype=np.float64)
    dshpg1 = np.empty(10, dtype=np.float64)
    dshpg2 = np.empty(10, dtype=np.float64)

    # dmatc = hooke(0, materialbyElement, 0, np.array(3 * [1.0]), np.array(6 * [1.0]), 0.0, 0.0, 0.0,
    #               210000.0)  # elastic material stiffness matrix

    # E = materialbyElement[0][0]  # Young's Modulus
    # Es = 210000.0
    # nu = materialbyElement[0][1]  # Poisson's Ratio
    # nu = 0.0  # for simplified concrete analysis
    # G = E / 2.0 / (1 + nu)  # shear modulus

    first = True

    for el, nodes in enumerate(elNodes):
        for i, nd in enumerate(nodes):
            co = nocoord[nd - 1]
            xlv0[i] = co[0]
            xlv1[i] = co[1]
            xlv2[i] = co[2]
            # xlv[i] = co

        elpos = 24 * el  # 4 integration points with 6 stress components each
        elpos_s = 12 * el  # 4 integration points with 3 stress components each
        elv = np.zeros(30, dtype=np.float64)  # element load vector, 10 nodes with 3 load components each

        for index, nd in enumerate(nodes):
            n3 = 3 * (nd - 1)
            i3 = 3 * index
            u10[i3] = du[n3]
            u10[i3 + 1] = du[n3 + 1]
            u10[i3 + 2] = du[n3 + 2]
        for i in range(4):

            sy = sig_yield[4 * el + i]
            ip = gp10[i]

            ipos1 = 4 * el + i
            ipos2 = 24 * el + 6 * i

            # print("Gauss point: ", 4 * el + i)

            # ------------------------------------------------------------------------------------------------------------
            xi = ip[0]
            et = ip[1]
            ze = ip[2]

            # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
            dshp0[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp0[1] = 4.0 * xi - 1.0
            dshp0[2] = 0.0
            dshp0[3] = 0.0
            dshp0[4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
            dshp0[5] = 4.0 * et
            dshp0[6] = -4.0 * et
            dshp0[7] = -4.0 * ze
            dshp0[8] = 4.0 * ze
            dshp0[9] = 0.0

            # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
            dshp1[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp1[1] = 0.0
            dshp1[2] = 4.0 * et - 1.0
            dshp1[3] = 0.0
            dshp1[4] = -4.0 * xi
            dshp1[5] = 4.0 * xi
            dshp1[6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
            dshp1[7] = -4.0 * ze
            dshp1[8] = 0.0
            dshp1[9] = 4.0 * ze

            # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
            dshp2[0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
            dshp2[1] = 0.0
            dshp2[2] = 0.0
            dshp2[3] = 4.0 * ze - 1.0
            dshp2[4] = -4.0 * xi
            dshp2[5] = 0.0
            dshp2[6] = -4.0 * et
            dshp2[7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
            dshp2[8] = 4.0 * xi
            dshp2[9] = 4.0 * et

            # local derivative of the global coordinates
            xs00 = xs01 = xs02 = xs10 = xs11 = xs12 = xs20 = xs21 = xs22 = 0.0
            for kb in range(10):
                xlv0kb = xlv0[kb]
                xlv1kb = xlv1[kb]
                xlv2kb = xlv2[kb]
                dshp0kb = dshp0[kb]
                dshp1kb = dshp1[kb]
                dshp2kb = dshp2[kb]
                xs00 += xlv0kb * dshp0kb
                xs01 += xlv0kb * dshp1kb
                xs02 += xlv0kb * dshp2kb
                xs10 += xlv1kb * dshp0kb
                xs11 += xlv1kb * dshp1kb
                xs12 += xlv1kb * dshp2kb
                xs20 += xlv2kb * dshp0kb
                xs21 += xlv2kb * dshp1kb
                xs22 += xlv2kb * dshp2kb

            # Jacobian
            xsj = (xs00 * xs11 * xs22 -
                   xs00 * xs12 * xs21 +
                   xs02 * xs10 * xs21 -
                   xs02 * xs11 * xs20 +
                   xs01 * xs12 * xs20 -
                   xs01 * xs10 * xs22)

            # global derivative of the local coordinates
            xsi00 = (xs11 * xs22 - xs21 * xs12) / xsj
            xsi01 = (xs02 * xs21 - xs01 * xs22) / xsj
            xsi02 = (xs01 * xs12 - xs02 * xs11) / xsj
            xsi10 = (xs12 * xs20 - xs10 * xs22) / xsj
            xsi11 = (xs00 * xs22 - xs02 * xs20) / xsj
            xsi12 = (xs10 * xs02 - xs00 * xs12) / xsj
            xsi20 = (xs10 * xs21 - xs20 * xs11) / xsj
            xsi21 = (xs20 * xs01 - xs00 * xs21) / xsj
            xsi22 = (xs00 * xs11 - xs10 * xs01) / xsj

            # global derivatives of the shape functions
            for jb in range(10):
                dshpg0[jb] = xsi00 * dshp0[jb] + xsi10 * dshp1[jb] + xsi20 * dshp2[jb]
                dshpg1[jb] = xsi01 * dshp0[jb] + xsi11 * dshp1[jb] + xsi21 * dshp2[jb]
                dshpg2[jb] = xsi02 * dshp0[jb] + xsi12 * dshp1[jb] + xsi22 * dshp2[jb]

            # strain interpolation matrix bm1 .. bm8
            for ib in range(10):
                d00 = dshpg0[ib]
                d10 = dshpg1[ib]
                d20 = dshpg2[ib]
                bm0[ib] = d00
                bm1[ib] = d10
                bm2[ib] = d20
                bm3[ib] = d10
                bm4[ib] = d00
                bm5[ib] = d20
                bm6[ib] = d00
                bm7[ib] = d20
                bm8[ib] = d10

            # ------------------------------------------------------------------------------------------

            # strain
            eps = np.zeros(6, dtype=np.float64)
            for j in range(10):
                j3 = 3 * j
                eps[0] += bm0[j] * u10[j3]
                eps[1] += bm1[j] * u10[j3 + 1]
                eps[2] += bm2[j] * u10[j3 + 2]
                eps[3] += bm3[j] * u10[j3] + bm4[j] * u10[j3 + 1]
                eps[4] += bm5[j] * u10[j3] + bm6[j] * u10[j3 + 2]
                eps[5] += bm7[j] * u10[j3 + 1] + bm8[j] * u10[j3 + 2]

            # print("eps: ", eps)

            dmatc = np.zeros((6, 6), dtype=np.float64)
            Ec = materialbyElement[el][0]
            dmatc[0][0] = dmatc[1][1] = dmatc[2][2] = Ec
            dmatc[3][3] = dmatc[4][4] = dmatc[5][5] = Ec / 2.0

            dmats = np.zeros((6, 6), dtype=np.float64)
            Es = 210000.0
            dmats[0, 0] = rhox * Es
            dmats[1, 1] = rhoy * Es
            dmats[2, 2] = rhoz * Es

            if model == 1:
                # elastic test stress
                for j in range(6):
                    tmp = sig[ipos2 + j]
                    for k in range(6):
                        tmp += (dmatc[j, k] + dmats[j, k]) * eps[k]
                    sig_test[j] = tmp
            elif model == 2:
                # elastic test stress concrete
                for j in range(6):
                    tmp = sigc[ipos2 + j]
                    for k in range(6):
                        tmp += dmatc[j, k] * eps[k]
                    sig_test[j] = tmp
                # elastic test stress steel
                for j in range(6):
                    tmp = sigs[ipos2 + j]
                    for k in range(6):
                        tmp += dmats[j, k] * eps[k]
                    sigs_test[j] = tmp

            # print("sig_test_concrete: ", sig_test_concrete)

            # sig_test_global[ipos2:ipos2 + 6] = sig_test

            prn = False

            sxx, syy, szz, sxy, syz, szx, psdir, act_c, act_s, depscc, depsct, depss = concrete(
                sig_test, rhox, rhoy, rhoz, fcd, sy, 0.0, step, iterat, prn, Ec, model)

            if model == 2:
                sigc_update[ipos2:ipos2 + 6] = sxx, syy, szz, sxy, syz, szx
                # ignore steel in compression
                sxxs = min(rhox * sy, max(sigs_test[0], 0.0))
                syys = min(rhoy * sy, max(sigs_test[1], 0.0))
                szzs = min(rhoz * sy, max(sigs_test[2], 0.0))
                sxx += sxxs
                syy += syys
                szz += szzs
                sigs_update[ipos2:ipos2 + 6] = sxxs, syys, szzs, 0.0, 0.0, 0.0
                dlx = max((sigs_test[0] - sxxs) / dmats[0, 0], 0.0)
                dly = max((sigs_test[1] - syys) / dmats[1, 1], 0.0)
                dlz = max((sigs_test[2] - szzs) / dmats[2, 2], 0.0)
                depss = max(dlx, dly, dlz)

            sig_update[ipos2:ipos2 + 6] = sxx, syy, szz, sxy, syz, szx

            # print(depsc)

            epscc_update[ipos1] = epscc[ipos1] + depscc
            epsct_update[ipos1] = epsct[ipos1] + depsct
            epss_update[ipos1] = epss[ipos1] + depss

            if act_c[0] == 1:
                ev[ipos2:ipos2 + 3] = psdir[:, 0]
                if act_c[1] == 1:
                    ev[ipos2 + 3:ipos2 + 6] = psdir[:, 1]
                else:
                    ev[ipos2 + 3:ipos2 + 6] = psdir[:, 2]
            else:
                if act_c[1] == 1:
                    ev[ipos2:ipos2 + 3] = psdir[:, 1]
                    ev[ipos2 + 3:ipos2 + 6] = psdir[:, 2]
                else:
                    ev[ipos2:ipos2 + 3] = psdir[:, 2]

            a_c[4 * el + i] = np.sum(act_c)
            a_s[12 * el + 3 * i:12 * el + 3 * i + 3] = act_s

            # calculate element load vectors
            ipxsj = ip[3] * abs(xsj)
            for j in range(10):
                j3 = 3 * j
                elv[j3] += (bm0[j] * sxx + bm3[j] * sxy + bm5[j] * syz) * ipxsj
                elv[j3 + 1] += (bm1[j] * syy + bm4[j] * sxy + bm7[j] * szx) * ipxsj
                elv[j3 + 2] += (bm2[j] * szz + bm6[j] * syz + bm8[j] * szx) * ipxsj

        # add element load vectors to the global load vector
        for i in range(10):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range(3):
                qin[iglob3 + k] += elv[i3 + k]

    return


@jit("types.UniTuple(float64, 6)(float64[::1], float64, float64, float64)", nopython=True, cache=True)
def vmises_original_optimised(sig_test, sig_yield, H, G):
    st0, st1, st2, st3, st4, st5 = sig_test
    p = (st0 + st1 + st2) / 3.0
    st0 -= p
    st1 -= p
    st2 -= p
    sig_mises = np.sqrt(1.5 * (st0 ** 2 + st1 ** 2 + st2 ** 2) +
                        3.0 * (st3 ** 2 + st4 ** 2 + st5 ** 2))

    if sig_yield > sig_mises:
        fac = 1.0
    else:
        fac = (1.0 - (1.0 - sig_yield / sig_mises) * 3.0 * G / (H + 3 * G))

    st0 = fac * st0 + p
    st1 = fac * st1 + p
    st2 = fac * st2 + p
    st3 = fac * st3
    st4 = fac * st4
    st5 = fac * st5

    sig_update = st0, st1, st2, st3, st4, st5

    return sig_update


@jit(nopython=True, cache=True, nogil=True)
def concrete(sig_test, rhox, rhoy, rhoz, fcd, sig_yield, nu, step, iterat, prn, E, model):
    a_c = np.array(3 * [0])
    a_s = np.array(3 * [1.0])

    depss = 0.0
    ps1, ps2, ps3, psdir = calculate_principal_stress_ip(sig_test)

    if model == 1:
        fy1 = (rhox * psdir[0, 0] ** 2 + rhoy * psdir[1, 0] ** 2 + rhoz * psdir[2, 0] ** 2) * sig_yield
        fy2 = (rhox * psdir[0, 1] ** 2 + rhoy * psdir[1, 1] ** 2 + rhoz * psdir[2, 1] ** 2) * sig_yield
        fy3 = (rhox * psdir[0, 2] ** 2 + rhoy * psdir[1, 2] ** 2 + rhoz * psdir[2, 2] ** 2) * sig_yield

        # ignore steel in compression and concrete in tension
        psr1 = min(fy1, max(ps1, -fcd))
        psr2 = min(fy2, max(ps2, -fcd))
        psr3 = min(fy3, max(ps3, -fcd))

        # plastic compressive strain increments concrete
        dl1 = -min((ps1 - psr1) / E, 0.0)
        dl2 = -min((ps2 - psr2) / E, 0.0)
        dl3 = -min((ps3 - psr3) / E, 0.0)
        depscc = max(dl1, dl2, dl3)

        # plastic tensile strain increments concrete
        dl1 = max((ps1 - psr1) / E, 0.0)
        dl2 = max((ps2 - psr2) / E, 0.0)
        dl3 = max((ps3 - psr3) / E, 0.0)
        depsct = max(dl1, dl2, dl3)
        # print("depsct: ", depsct)

        # plastic tensile strain increments steel
        dlx = dl1 * psdir[0, 0] ** 2 + dl2 * psdir[0, 1] ** 2 + dl3 * psdir[0, 2] ** 2
        dly = dl1 * psdir[1, 0] ** 2 + dl2 * psdir[1, 1] ** 2 + dl3 * psdir[1, 2] ** 2
        dlz = dl1 * psdir[2, 0] ** 2 + dl2 * psdir[2, 1] ** 2 + dl3 * psdir[2, 2] ** 2
        depss = max(dlx, dly, dlz)
        # print("depss: ", depss)
        # print("dlx, dly, dlz: ", dlx, dly, dlz)


    elif model == 2:
        # concrete stress
        psr1 = min(0.0, max(ps1, -fcd))
        psr2 = min(0.0, max(ps2, -fcd))
        psr3 = min(0.0, max(ps3, -fcd))

        # plastic compressive strain increments concrete
        dl1 = -min((ps1 - psr1) / E, 0.0)
        dl2 = -min((ps2 - psr2) / E, 0.0)
        dl3 = -min((ps3 - psr3) / E, 0.0)

        # plastic tensile strain increments concrete
        dl1 = max((ps1 - psr1) / E, 0.0)
        dl2 = max((ps2 - psr2) / E, 0.0)
        dl3 = max((ps3 - psr3) / E, 0.0)
        depsct = max(dl1, dl2, dl3)

        if ps1 > 1.0e-3: a_c[0] = 1
        if ps2 > 1.0e-3: a_c[1] = 1
        if ps3 > 1.0e-3: a_c[2] = 1
        # print(ps1,ps2,ps3)

        # equivalent plastic compressive strain increment concrete
        # depsc = np.sqrt(2 / 3 * (dl1 ** 2 + dl2 ** 2 + dl3 ** 2))
        depsc = max(dl1, dl2, dl3)

        depss = 0.0

    sigxx_c, sigyy_c, sigzz_c, sigxy, sigyz, sigzx = cartesian(psr1, psr2, psr3, psdir)

    return sigxx_c, sigyy_c, sigzz_c, sigxy, sigyz, sigzx, psdir, a_c, a_s, depscc, depsct, depss


@jit(nopython=True, cache=True, nogil=True)
def rankine_modified(ps1, ps2, ps3, nu, step, iterat):
    rankine_elastic = True
    psr = np.array(3 * [0.])
    ps = np.array([ps1, ps2, ps3])

    f1e = ps[0]
    f2e = ps[1]
    f3e = ps[2]

    if max(ps) > 0.0:
        psr[0] = 0.0
        psr[1] = ps[1] - nu / (1.0 - nu) * ps[0]
        psr[2] = ps[2] - nu / (1.0 - nu) * ps[0]
        # if step == 11 and iterat == 0:
        # print(">>>> REGULAR")
        # print(">>>> psr: ", psr[0], psr[1], psr[2])

        if psr[1] > 1.0e-6:
            dl1 = ((1.0 - nu) * f1e - nu * f2e)
            dl2 = (-nu * f1e + (1.0 - nu) * f2e)
            psr[0] = ps[0] - (dl1 * (1.0 - nu) + dl2 * nu) / (1.0 - 2.0 * nu)
            psr[1] = ps[1] - (dl1 * nu + dl2 * (1.0 - nu)) / (1.0 - 2.0 * nu)
            psr[2] = ps[2] - (dl1 * nu + dl2 * nu) / (1.0 - 2.0 * nu)
            # if step == 11 and iterat == 0:
            # print(">>>> RIDGE12")
            # print(">>>> psr: ", psr[0], psr[1], psr[2])

        if psr[2] > 1.0e-6:
            dl1 = ((1.0 - nu) * f1e - nu * f3e)
            dl3 = (-nu * f1e + (1.0 - nu) * f3e)
            psr[0] = ps[0] - (dl1 * (1.0 - nu) + dl3 * nu) / (1.0 - 2.0 * nu)
            psr[1] = ps[1] - (dl1 * nu + dl3 * nu) / (1.0 - 2.0 * nu)
            psr[2] = ps[2] - (dl1 * nu + dl3 * (1.0 - nu)) / (1.0 - 2.0 * nu)
            # if step == 11 and iterat == 0:
            # print(">>>> RIDGE23")
            # print(">>>> psr: ", psr[0], psr[1], psr[2])

        if psr[1] > 1.0e-6 or psr[2] > 1.0e-6:
            # if step == 11 and iterat == 0: print(">>>> APEX")
            psr[0] = 0.0
            psr[1] = 0.0
            psr[2] = 0.0

        ps1 = psr[0]
        ps2 = psr[1]
        ps3 = psr[2]

    if step == 0 and iterat == 0:
        print("stresses at the end of rankine: ", ps1, ps2, ps3)
    return ps1, ps2, ps3


@jit(nopython=True, cache=True, nogil=True)
def rankine(ps1, ps2, ps3, ts1, ts2, ts3, nu, first):
    rankine_elastic = True
    psr = np.array(3 * [0.])
    ps = np.array([ps1, ps2, ps3])
    if first:
        print(">>>> ps: ", ps, "ts: ", ts, "ps: ", ps)
        # print("ts: ", ts)
        # print("ps: ", ps)

    if max(ps) > 0.0:
        rankine_elastic = False
        si = np.argsort(ps)[::-1]
        psr[si[0]] = ps[si[0]] - ps[si[0]]
        psr[si[1]] = ps[si[1]] - nu / (1.0 - nu) * ps[si[0]]
        psr[si[2]] = ps[si[2]] - nu / (1.0 - nu) * ps[si[0]]
        # print("si[0]: ", si[0])
        # print("psr[si[0]]: ", psr[si[0]])
        # print("psr[si[1]]: ", psr[si[1]])
        # print("psr[si[2]]: ", psr[si[2]])

        f2r = psr[si[1]] - ts[si[1]]
        f3r = psr[si[2]] - ts[si[2]]

        # print("f2r: ", f2r)
        # print("f3r: ", f3r)

        f1e = ps[si[0]] - ts[si[0]]
        f2e = ps[si[1]] - ts[si[1]]
        f3e = ps[si[2]] - ts[si[2]]

        if first:
            print(">>>> f2r: ", f2r)
            print(">>>> f3r: ", f3r)

        apex = (f2r > 0.0 and f3r > 0.0)  # test for apex

        if apex:
            ps1 = ts[si[0]]
            ps2 = ts[si[1]]
            ps3 = ts[si[2]]
            if first: print(">>>> APEX")
        else:
            ridge12 = f2r > 0.0  # test for ps1 = ps2 ridge
            ridge23 = f3r > 0.0  # test for ps2 = ps3 ridge

            # print("ridge12: ", ridge12)
            # print("ridge23: ", ridge23)

            if ridge12:
                if first: print("RIDGE12")
                dl1 = ((1.0 - nu) * f1e - nu * f2e)
                dl2 = (-nu * f1e + (1.0 - nu) * f2e)
                ps1 = ps[si[0]] - (dl1 * (1.0 - nu) + dl2 * nu) / (1.0 - 2.0 * nu)
                ps2 = ps[si[1]] - (dl1 * nu + dl2 * (1.0 - nu)) / (1.0 - 2.0 * nu)
                ps3 = ps[si[2]] - (dl1 * nu + dl2 * nu) / (1.0 - 2.0 * nu)
            elif ridge23:
                if first: print("RIDGE23")
                dl1 = ((1.0 - nu) * f1e - nu * f3e)
                dl3 = (-nu * f1e + (1.0 - nu) * f3e)
                ps1 = ps[si[0]] - (dl1 * (1.0 - nu) + dl3 * nu) / (1.0 - 2.0 * nu)
                ps2 = ps[si[1]] - (dl1 * nu + dl3 * nu) / (1.0 - 2.0 * nu)
                ps3 = ps[si[2]] - (dl1 * nu + dl3 * (1.0 - nu)) / (1.0 - 2.0 * nu)
            else:
                if first:
                    print("REGULAR")
                    print("si: ", si)

                ps1 = psr[0]
                ps2 = psr[1]
                ps3 = psr[2]
                # ps1 = psr[si[0]]
                # ps2 = psr[si[1]]
                # ps3 = psr[si[2]]

    # print("stresses at the end of rankine: ", ps1, ps2, ps3)
    return ps1, ps2, ps3


@jit(nopython=True, cache=True, nogil=True)
def tresca(ps1, ps2, ps3, cf):
    femc = (ps1 - ps3) - 2.0 * cf

    if femc > 0.0:
        ps1r = ps1 - femc / 2.0
        ps2r = ps2
        ps3r = ps3 + femc / 2.0

        ridge12 = ps2r > 1.0000001 * ps1r  # test for ps1 = ps2 ridge
        ridge23 = 1.0000001 * ps2r < ps3r  # test for ps2 = ps3 ridge

        if ridge12:  # triaxial compression
            ridge = 1
            femctc = ((ps2 - ps3) - 2.0 * cf)
            ps1 -= (2.0 * femc - femctc) / 3.0
            ps2 -= (2.0 * femctc - femc) / 3.0
            ps3 += (femc + femctc) / 3.0

        elif ridge23:  # triaxial extension
            ridge = 2
            femcte = ((ps1 - ps2) - 2.0 * cf)
            ps1 -= (femc + femcte) / 3.0
            ps2 += (2.0 * femcte - femc) / 3.0
            ps3 += (2.0 * femc - femcte) / 3.0

        else:  # regular return
            ridge = 0
            ps1 = ps1r
            ps2 = ps2r
            ps3 = ps3r

    return ps1, ps2, ps3


@jit(nopython=True, cache=True, nogil=True)
def calculate_principal_stress_ip(stress_tensor):
    s11 = stress_tensor[0]  # Sxx
    s22 = stress_tensor[1]  # Syy
    s33 = stress_tensor[2]  # Szz
    s12 = stress_tensor[3]  # Sxy
    s31 = stress_tensor[4]  # Szx
    s23 = stress_tensor[5]  # Syz
    sigma = np.array([
        [s11, s12, s31],
        [s12, s22, s23],
        [s31, s23, s33]
    ])

    eigenvalues, eigenvectors = np.linalg.eigh(sigma)

    eigenvalues = eigenvalues
    eigenvectors = eigenvectors

    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    return (eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvectors)


@jit(nopython=True, cache=True, nogil=True)
def cartesian(ps1, ps2, ps3, psdir):
    sig_update = np.zeros(6, dtype=types.float64)  # updated stresses
    sig_update[0] = psdir[0, 0] ** 2 * ps1 + psdir[0, 1] ** 2 * ps2 + psdir[0, 2] ** 2 * ps3
    sig_update[1] = psdir[1, 0] ** 2 * ps1 + psdir[1, 1] ** 2 * ps2 + psdir[1, 2] ** 2 * ps3
    sig_update[2] = psdir[2, 0] ** 2 * ps1 + psdir[2, 1] ** 2 * ps2 + psdir[2, 2] ** 2 * ps3
    sig_update[3] = psdir[0, 0] * psdir[1, 0] * ps1 + psdir[0, 1] * psdir[1, 1] * ps2 + psdir[0, 2] * \
                    psdir[1, 2] * ps3
    sig_update[4] = psdir[0, 0] * psdir[2, 0] * ps1 + psdir[0, 1] * psdir[2, 1] * ps2 + psdir[0, 2] * \
                    psdir[2, 2] * ps3
    sig_update[5] = psdir[1, 0] * psdir[2, 0] * ps1 + psdir[1, 1] * psdir[2, 1] * ps2 + psdir[1, 2] * \
                    psdir[2, 2] * ps3

    return sig_update


# map stresses to nodes
@jit(nopython=True, cache=True)
def mapStresses(averaged, elNodes, nocoord, sig, epscc, epsct, sigvm, epss, noce, sig_yield, pdfcd):
    # map maps corner node stresses to all tet10 nodes

    map_inter = np.array([[0.5, 0.5, 0.0, 0.0],
                          [0.0, 0.5, 0.5, 0.0],
                          [0.5, 0.0, 0.5, 0.0],
                          [0.5, 0.0, 0.0, 0.5],
                          [0.0, 0.5, 0.0, 0.5],
                          [0.0, 0.0, 0.5, 0.5]])

    tet10stress = np.zeros((len(nocoord), 6), dtype=np.float64)
    tet10epss = np.zeros(len(nocoord), dtype=np.float64)
    tet10epscc = np.zeros(len(nocoord), dtype=np.float64)
    tet10epsct = np.zeros(len(nocoord), dtype=np.float64)
    tet10svm = np.zeros(len(nocoord), dtype=np.float64)
    tet10pdfcd = np.zeros(len(nocoord), dtype=np.float64)

    tet4stress = sig.reshape((len(elNodes), 4, 6))
    tet4epss = epss.reshape((len(elNodes), 4))
    tet4epscc = epscc.reshape((len(elNodes), 4))
    tet4epsct = epsct.reshape((len(elNodes), 4))
    tet4svm = sigvm.reshape((len(elNodes), 4))
    tet4pdfcd = pdfcd.reshape((len(elNodes), 4))

    # averaged nodal stresses
    for el, nodes in enumerate(elNodes):
        corners = nodes[0:4]
        tet10stress[corners - 1] += tet4stress[el] / noce[corners - 1].reshape((4, 1))

    if averaged:
        # averaged nodal scalars
        for el, nodes in enumerate(elNodes):
            corners = nodes[0:4]
            tet10epss[corners - 1] += tet4epss[el] / noce[corners - 1]
            tet10epscc[corners - 1] += tet4epscc[el] / noce[corners - 1]
            tet10epsct[corners - 1] += tet4epsct[el] / noce[corners - 1]
            tet10svm[corners - 1] += tet4svm[el] / noce[corners - 1]
            tet10pdfcd[corners - 1] += tet4pdfcd[el] / noce[corners - 1]
    else:
        # unaveraged nodal scalars
        for el, nodes in enumerate(elNodes):
            corners = nodes[0:4]
            tet10epss[corners - 1] = np.fmax(tet10epss[corners - 1], tet4epss[el])
            tet10epscc[corners - 1] = np.fmax(tet10epscc[corners - 1], tet4epscc[el])
            tet10epsct[corners - 1] = np.fmax(tet10epsct[corners - 1], tet4epsct[el])
            tet10svm[corners - 1] = np.fmax(tet10svm[corners - 1], tet4svm[el])
            tet10pdfcd[corners - 1] = np.fmin(tet10pdfcd[corners - 1], tet4pdfcd[el])

    # results intermediate nodes
    for el, nodes in enumerate(elNodes):
        nd_corner = nodes[0:4]
        nd_inter = nodes[4:10]
        tet10stress[nd_inter - 1] = np.dot(map_inter, tet10stress[nd_corner - 1])
        tet10epss[nd_inter - 1] = np.dot(map_inter, tet10epss[nd_corner - 1])
        tet10epscc[nd_inter - 1] = np.dot(map_inter, tet10epscc[nd_corner - 1])
        tet10epsct[nd_inter - 1] = np.dot(map_inter, tet10epsct[nd_corner - 1])
        tet10svm[nd_inter - 1] = np.dot(map_inter, tet10svm[nd_corner - 1])
        tet10pdfcd[nd_inter - 1] = np.dot(map_inter, tet10pdfcd[nd_corner - 1])

    return tet10stress, tet10epss, tet10epscc, tet10epsct, tet10svm, tet10pdfcd


# fill resultobject with results
def pasteResults(doc, elNodes, nocoord, dis, tet10stress, tet10epss, tet10epscc, tet10rho):
    return_code = 0
    analysis = doc.getObject("Analysis")

    nn = len(nocoord)  # number of nodes

    resVol = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]

    # VOLUME MESH START
    elements_tetra10 = {}
    mode_results_vol = {}

    nodes = range(1, nn + 1)
    coordinates = map(App.Vector, nocoord)
    displacements = map(App.Vector, np.array_split(dis, nn))
    volnodes = dict(zip(nodes, coordinates))
    mode_disp_vol = dict(zip(nodes, displacements))

    for index, elem in enumerate(elNodes):
        elements_tetra10[index + 1] = (
            elem[0], elem[2], elem[1], elem[3], elem[6], elem[5], elem[4], elem[7], elem[9], elem[8])

    mode_results_vol['disp'] = mode_disp_vol

    results = [mode_results_vol]

    mvol = {
        'Nodes': volnodes,
        'Seg2Elem': {},
        'Seg3Elem': {},
        'Tria3Elem': {},
        'Tria6Elem': {},
        'Quad4Elem': {},
        'Quad8Elem': {},
        'Tetra4Elem': {},
        'Tetra10Elem': elements_tetra10,
        'Hexa8Elem': {},
        'Hexa20Elem': {},
        'Penta6Elem': {},
        'Penta15Elem': {},
        'Results': results
    }

    meshvol = itf.make_femmesh(mvol)

    result_mesh_object_1 = ObjectsFem.makeMeshResult(doc, 'Result_Mesh_Volume')
    result_mesh_object_1.FemMesh = meshvol

    resVol.DisplacementVectors = [App.Vector(dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]) for n in range(nn)]
    resVol.DisplacementLengths = [np.linalg.norm([dis[3 * n], dis[3 * n + 1], dis[3 * n + 2]]) for n in range(nn)]
    resVol.NodeStressXX = tet10stress.T[0].T.tolist()
    resVol.NodeStressYY = tet10stress.T[1].T.tolist()
    resVol.NodeStressZZ = tet10stress.T[2].T.tolist()
    resVol.NodeStressXY = tet10stress.T[3].T.tolist()
    resVol.NodeStressXZ = tet10stress.T[4].T.tolist()
    resVol.NodeStressYZ = tet10stress.T[5].T.tolist()

    try:
        resVol.CriticalStrainRatio = tet10epscc.T.tolist()  # FreeCAD 0.21.0 and higher - works for export to VTK
        resVol.Peeq = tet10epss.T.tolist()
        resVol.Temperature = tet10epscc.T.tolist()
        resVol.vonMises = (100 * tet10rho[:, 0].T).tolist()
        resVol.PrincipalMax = (100 * tet10rho[:, 1].T).tolist()
        resVol.PrincipalMin = (100 * tet10rho[:, 2]).T.tolist()

    except:
        resVol.Peeq = tet10epscc.T.tolist()  # a hack for FreeCAD 0.20.x - store the critical strain ratio in the PEEQ output
        resVol.Temperature = tet10epscc.T.tolist()
        resVol.vonMises = (100 * tet10rho[:, 0].T).tolist()
        resVol.PrincipalMax = (100 * tet10rho[:, 1].T).tolist()
        resVol.PrincipalMin = (100 * tet10rho[:, 2]).T.tolist()

    print("max(tet10rho[:, 0]): ", max(tet10rho[:, 0]))
    print("max(tet10rho[:, 1]): ", max(tet10rho[:, 1]))
    print("max(tet10rho[:, 2]): ", max(tet10rho[:, 2]))

    resVol.Stats = [min(dis[0::3]), max(dis[0::3]),
                    min(dis[1::3]), max(dis[1::3]),
                    min(dis[2::3]), max(dis[2::3]),
                    min(resVol.DisplacementLengths), max(resVol.DisplacementLengths),
                    min(resVol.vonMises), max(resVol.vonMises),
                    min(resVol.PrincipalMax), max(resVol.PrincipalMax),
                    0.0, 0.0,
                    min(resVol.PrincipalMin), max(resVol.PrincipalMin),
                    0.0, 0.0,
                    min(resVol.Peeq), max(resVol.Peeq),
                    min(resVol.Temperature), max(resVol.Temperature),
                    0.0, 0.0,
                    0.0, 0.0]

    resVol.Mesh = result_mesh_object_1
    resVol.NodeNumbers = [int(key) for key in resVol.Mesh.FemMesh.Nodes.keys()]

    resVol = itf.fill_femresult_mechanical(resVol, results)

    # VOLUME MESH FINISH

    # Add von Mises Stress and Plastic Strain Ratio to the results
    # rt.add_von_mises(resVol)
    # rt.add_principal_stress_std(resVol)
    #
    # rt.fill_femresult_stats(resVol)

    trm._TaskPanel.result_obj = resVol
    trm._TaskPanel.mesh_obj = meshvol

    # print(dir(trm._TaskPanel))

    doc.recompute()

    return resVol


def calcSum(Edge_Elements, Face_Elements, mesh, CSR, peeq, svm):
    gp10, gp6, gp2 = gaussPoints()

    coor = mesh.Nodes

    edge_length = []
    edge_peeq = []
    edge_CSR = []
    edge_svm = []

    for index, edge in enumerate(Edge_Elements):
        edge_length.append(0.0)
        edge_peeq.append(0.0)
        edge_CSR.append(0.0)
        edge_svm.append(0.0)
        for element in edge:
            xle = []
            for node in element:
                xle.append([coor[node].x, coor[node].y, coor[node].z])
            # integrate variable
            xle = np.array(xle).T
            for gp in gp2:
                xi = gp[0]
                xsj, shp = shape2lin(xi, xle)
                for i in range(3):
                    nd = element[i] - 1
                    dl = shp[i] * abs(xsj) * gp[1]
                    edge_length[index] += dl
                    edge_peeq[index] += peeq[nd] * dl
                    edge_CSR[index] += CSR[nd] * dl
                    edge_svm[index] += svm[nd] * dl
        Length = edge_length[index]
        if Length > 0.0:
            edge_peeq[index] /= Length
            edge_CSR[index] /= Length
            edge_svm[index] /= Length

    face_area = []
    face_peeq = []
    face_CSR = []
    face_svm = []

    for index, face in enumerate(Face_Elements):
        face_area.append(0.0)
        face_peeq.append(0.0)
        face_CSR.append(0.0)
        face_svm.append(0.0)
        for element in face:
            xlf = []
            for node in element:
                xlf.append([coor[node].x, coor[node].y, coor[node].z])
            # integrate variable
            xlf = np.array(xlf).T
            for gp in gp6:
                xi = gp[0]
                et = gp[1]
                xsj, shp, bmatS, xx, xt, xp = shape6tri(xi, et, xlf)
                for i in range(6):
                    nd = element[i] - 1
                    dA = shp[i] * abs(xsj) * gp[2]
                    face_area[index] += dA
                    face_peeq[index] += peeq[nd] * dA
                    face_CSR[index] += CSR[nd] * dA
                    face_svm[index] += svm[nd] * dA
        Area = face_area[index]
        if Area > 0.0:
            face_peeq[index] /= Area
            face_CSR[index] /= Area
            face_svm[index] /= Area

    return edge_length, edge_peeq, edge_CSR, edge_svm, face_area, face_peeq, face_CSR, face_svm


def exportVTK(elNodes, nocoord, dis, tet10stress, tet10epss, tet10epscc, tet10epsct, tet10svm, tet10pdfcd, fy, file):
    padding = np.full(len(elNodes), 10, dtype=int)
    cells = np.vstack((padding, (elNodes - 1).T)).T

    celltypes = np.full(len(elNodes), fill_value=CellType.QUADRATIC_TETRA, dtype=np.uint8)

    grid = pv.UnstructuredGrid(cells, celltypes, nocoord)
    grid.point_data["Concrete Plastic Compressive Strain\n"] = tet10epscc.flatten(order="F")
    grid.point_data["Concrete Plastic Tensile Strain\n"] = tet10epsct.flatten(order="F")
    grid.point_data["Steel Plastic Tensile Strain\n"] = tet10epss.flatten(order="F")
    # grid.point_data["von Mises Stress\n"] = tet10svm.flatten(order="F")
    grid.point_data["p/fcd\n"] = tet10pdfcd.flatten(order="F")

    displacement = dis.reshape((len(nocoord), 3))

    stress = tet10stress.reshape((len(nocoord), 6))

    grid.point_data['Displacement'] = displacement

    grid.point_data['Stress Tensor'] = stress

    tet10s1, tet10s2, tet10s3, sv1, sv2, sv3 = calculate_principal_stress(tet10stress)

    tet10rho = calculate_rho(tet10stress, fy)

    grid.point_data["Major Principal Stress\n"] = tet10s1.flatten(order="F")
    grid.point_data["Intermediate Principal Stress\n"] = tet10s2.flatten(order="F")
    grid.point_data["Minor Principal Stress\n"] = tet10s3.flatten(order="F")
    grid.point_data['Major Principal Stress Vector'] = sv1
    grid.point_data['Intermediate Principal Stress Vector'] = sv2
    grid.point_data['Minor Principal Stress Vector'] = sv3

    grid.point_data['Reinforcement Ratio x'] = tet10rho[:, 0].flatten(order="F")
    grid.point_data['Reinforcement Ratio y'] = tet10rho[:, 1].flatten(order="F")
    grid.point_data['Reinforcement Ratio z'] = tet10rho[:, 2].flatten(order="F")

    pv.save_meshio(file, grid)

    return tet10rho


@jit(nopython=True, cache=True, nogil=True)
def calculate_principal_stress(tet10stress):
    s1 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    s2 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    s3 = np.zeros((len(tet10stress), 3), dtype=np.float64)
    tet10s1 = np.zeros(len(tet10stress), dtype=np.float64)
    tet10s2 = np.zeros(len(tet10stress), dtype=np.float64)
    tet10s3 = np.zeros(len(tet10stress), dtype=np.float64)

    for index, stress in enumerate(tet10stress):
        s11 = stress[0]  # Sxx
        s22 = stress[1]  # Syy
        s33 = stress[2]  # Szz
        s12 = stress[3]  # Sxy
        s31 = stress[4]  # Szx
        s23 = stress[5]  # Syz
        sigma = np.array([
            [s11, s12, s31],
            [s12, s22, s23],
            [s31, s23, s33]
        ])

        # print(sigma)

        eigenvalues, eigenvectors = np.linalg.eigh(sigma)

        eigenvalues = eigenvalues.real
        eigenvectors = eigenvectors.real

        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        tet10s1[index] = eigenvalues[0]
        tet10s2[index] = eigenvalues[1]
        tet10s3[index] = eigenvalues[2]

        s1[index] = eigenvalues[0] * eigenvectors[:, 0]
        s2[index] = eigenvalues[1] * eigenvectors[:, 1]
        s3[index] = eigenvalues[2] * eigenvectors[:, 2]

    return tet10s1, tet10s2, tet10s3, s1, s2, s3


def calculate_rho(tet10stress, fy):
    #
    #   Calculation of Reinforcement Ratios and
    #   Concrete Stresses according to http://heronjournal.nl/53-4/3.pdf
    #           - See post:
    #             https://forum.freecadweb.org/viewtopic.php?f=18&t=28821
    #                   fy: factored yield strength of reinforcement bars
    #

    tet10rho = np.zeros((len(tet10stress), 3), dtype=np.float64)

    for index, stress_tensor in enumerate(tet10stress):
        rmin = 1.0e9
        eqmin = 14

        sxx = stress_tensor[0]
        syy = stress_tensor[1]
        szz = stress_tensor[2]
        sxy = stress_tensor[3]
        syz = stress_tensor[5]
        sxz = stress_tensor[4]

        rhox = np.zeros(15)
        rhoy = np.zeros(15)
        rhoz = np.zeros(15)

        #    i1=sxx+syy+szz NOT USED
        #    i2=sxx*syy+syy*szz+szz*sxx-sxy**2-sxz**2-syz**2 NOT USED
        i3 = (sxx * syy * szz + 2 * sxy * sxz * syz - sxx * syz ** 2
              - syy * sxz ** 2 - szz * sxy ** 2)

        #    Solution (5)
        d = (sxx * syy - sxy ** 2)
        if d != 0.:
            rhoz[0] = i3 / d / fy

        #    Solution (6)
        d = (sxx * szz - sxz ** 2)
        if d != 0.:
            rhoy[1] = i3 / d / fy

        #    Solution (7)
        d = (syy * szz - syz ** 2)
        if d != 0.:
            rhox[2] = i3 / d / fy

        #    Solution (9)
        if sxx != 0.:
            fc = sxz * sxy / sxx - syz
            fxy = sxy ** 2 / sxx
            fxz = sxz ** 2 / sxx

            #    Solution (9+)
            rhoy[3] = syy - fxy + fc
            rhoy[3] /= fy
            rhoz[3] = szz - fxz + fc
            rhoz[3] /= fy

            #    Solution (9-)
            rhoy[4] = syy - fxy - fc
            rhoy[4] /= fy
            rhoz[4] = szz - fxz - fc
            rhoz[4] /= fy

        #   Solution (10)
        if syy != 0.:
            fc = syz * sxy / syy - sxz
            fxy = sxy ** 2 / syy
            fyz = syz ** 2 / syy

            # Solution (10+)
            rhox[5] = sxx - fxy + fc
            rhox[5] /= fy
            rhoz[5] = szz - fyz + fc
            rhoz[5] /= fy

            # Solution (10-)vm
            rhox[6] = sxx - fxy - fc

            rhox[6] /= fy
            rhoz[6] = szz - fyz - fc
            rhoz[6] /= fy

        # Solution (11)
        if szz != 0.:
            fc = sxz * syz / szz - sxy
            fxz = sxz ** 2 / szz
            fyz = syz ** 2 / szz

            # Solution (11+)
            rhox[7] = sxx - fxz + fc
            rhox[7] /= fy
            rhoy[7] = syy - fyz + fc
            rhoy[7] /= fy

            # Solution (11-)
            rhox[8] = sxx - fxz - fc
            rhox[8] /= fy
            rhoy[8] = syy - fyz - fc
            rhoy[8] /= fy

        # Solution (13)
        rhox[9] = (sxx + sxy + sxz) / fy
        rhoy[9] = (syy + sxy + syz) / fy
        rhoz[9] = (szz + sxz + syz) / fy

        # Solution (14)
        rhox[10] = (sxx + sxy - sxz) / fy
        rhoy[10] = (syy + sxy - syz) / fy
        rhoz[10] = (szz - sxz - syz) / fy

        # Solution (15)
        rhox[11] = (sxx - sxy - sxz) / fy
        rhoy[11] = (syy - sxy + syz) / fy
        rhoz[11] = (szz - sxz + syz) / fy

        # Solution (16)
        rhox[12] = (sxx - sxy + sxz) / fy
        rhoy[12] = (syy - sxy - syz) / fy
        rhoz[12] = (szz + sxz - syz) / fy

        # Solution (17)
        if syz != 0.:
            rhox[13] = (sxx - sxy * sxz / syz) / fy
        if sxz != 0.:
            rhoy[13] = (syy - sxy * syz / sxz) / fy
        if sxy != 0.:
            rhoz[13] = (szz - sxz * syz / sxy) / fy

        for ir in range(0, rhox.size):

            if rhox[ir] >= -1.e-10 and rhoy[ir] >= -1.e-10 and rhoz[ir] > -1.e-10:

                # Concrete Stresses
                scxx = sxx - rhox[ir] * fy
                scyy = syy - rhoy[ir] * fy
                sczz = szz - rhoz[ir] * fy
                ic1 = (scxx + scyy + sczz)
                ic2 = (scxx * scyy + scyy * sczz + sczz * scxx - sxy ** 2
                       - sxz ** 2 - syz ** 2)
                ic3 = (scxx * scyy * sczz + 2 * sxy * sxz * syz - scxx * syz ** 2
                       - scyy * sxz ** 2 - sczz * sxy ** 2)

                if ic1 <= 1.e-6 and ic2 >= -1.e-6 and ic3 <= 1.0e-6:

                    rsum = rhox[ir] + rhoy[ir] + rhoz[ir]

                    if rsum < rmin and rsum > 0.:
                        rmin = rsum
                        eqmin = ir

        tet10rho[index] = [rhox[eqmin], rhoy[eqmin], rhoz[eqmin]]

    return tet10rho


def calculate_mohr_coulomb(prin1, prin3, phi, fck):
    #
    #             Calculation of Mohr Coulomb yield criterion to judge
    #             concrete curshing and shear failure
    #                   phi: angle of internal friction
    #                   fck: factored compressive strength of the matrix material (usually concrete)
    #

    coh = fck * (1 - np.sin(phi)) / 2 / np.cos(phi)

    mc_stress = ((prin1 - prin3) + (prin1 + prin3) * np.sin(phi)
                 - 2. * coh * np.cos(phi))

    if mc_stress < 0.:
        mc_stress = 0.

    return mc_stress
