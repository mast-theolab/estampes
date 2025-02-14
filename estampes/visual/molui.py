"""Build user interfaces for molecular visualization.

Provides classes to build simple or more advanced interfaces to
visualize molecules and associated properties.
"""
import os
import typing as tp

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base.types import Type1Vib, TypeAtCrdM, TypeAtLabM, TypeBondsM, \
    TypeColor
from estampes.visual.molview import Molecule
from estampes.visual.vibview import VibMode
from estampes.visual.povrender import POVBuilder, MSG_POV_CMD_HTML, MSG_POV_WIN


class MolWin(Qt3DExtras.Qt3DWindow):
    """The Molecule Window class (Qt3D).

    Generates a Qt3D window instance for the visualization of
    molecule(s).
    """

    click_mol = QtCore.Signal(list)

    def __init__(self, nmols: int,
                 atlabs: TypeAtLabM,
                 atcrds: TypeAtCrdM,
                 bonds: TypeBondsM,
                 model: tp.Optional[str] = None,
                 material: tp.Optional[str] = None,
                 col_bond_as_atom: bool = False,
                 molcols: tp.Optional[TypeColor] = None,
                 **extras):
        """Build an instance of MolWin.

        Builds an instance of MolWin.

        ----------
        nmols
            Number of molecules stored in `atlabs`, `atcrds` and `bonds`.
        atlabs
            Atomic labels.
            If `nmols>1`, list of lists.
        atcrds
            3-tuples with atomic coordinates, in Ang.
            If `nmols>1`, list of lists.
        bonds
            2-tuples listing connected atoms.
            If `nmols>1`, list of lists.
        model
            Representation model of the molecule.
        material
            Material simulated for the molecule.
        col_bond_as_atom
            If true, bonds are colored based on the connected atoms
        molcols
            If not `None`, color of the each molecule.
        extras
            Extra keywords.  Supported:

            `fname`: filename where data were stored.
            `fnames`: sequence of filenames where data were stored.
        """
        super(MolWin, self).__init__()

        self.__vib = None
        self.__mol = None
        self.__nmols = 0

        if 'fnames' in extras:
            fname = extras['fnames']
        elif 'fname' in extras:
            fname = extras['fname']
        else:
            fname = None
        if isinstance(fname, str):
            if nmols > 1:
                raise ValueError('Expected list of filenames')
            self.__fnames = fname
        elif isinstance(fname, (tuple, list)):
            if nmols == 1:
                self.__fnames = fname[0]
            else:
                self.__fnames = fname[:]
        else:
            self.__fnames = None

        # Camera
        self.__cam = self.camera()
        self.__cam.lens().setPerspectiveProjection(45, 16 / 9, 0.1, 1000)
        self.__cam.setPosition(QtGui.QVector3D(0, 0, 40))
        self.__cam.setUpVector(QtGui.QVector3D(0, 1, 0))
        self.__cam.setViewCenter(QtGui.QVector3D(0, 0, 0))

        # For camera controls
        self.rootEntity = Qt3DCore.QEntity()
        Molecule.reset_counter()
        if nmols == 1:
            self.__nmols = 1
            self.__mol = Molecule(atlabs, atcrds, bonds, model, material,
                                  col_bond_as_atom, molcols, self.rootEntity)
            self.__mol.addMouse(self.camera)
            self.__mol.click_molatom.connect(self.atom_clicked)
        else:
            self.__nmols = nmols
            self.__mol = []
            for i in range(nmols):
                if molcols is None:
                    molcol = None
                else:
                    molcol = molcols[i]
                self.__mol.append(Molecule(atlabs[i], atcrds[i], bonds[i],
                                           model, material, col_bond_as_atom,
                                           molcol, self.rootEntity))
                self.__mol[-1].addMouse(self.__cam)
                self.__mol[-1].click_molatom.connect(self.atom_clicked)
        self.camController = Qt3DExtras.QOrbitCameraController(self.rootEntity)
        self.camController.setLinearSpeed(50)
        self.camController.setLookSpeed(180)
        self.camController.setCamera(self.__cam)
        self.obj_light = Qt3DCore.QEntity(self.__cam)
        self.camLight = Qt3DRender.QPointLight(self.obj_light)
        self.cam_tvec = Qt3DCore.QTransform()
        self.cam_tvec.setTranslation(QtGui.QVector3D(0, 50, 100))
        self.obj_light.addComponent(self.camLight)
        self.obj_light.addComponent(self.cam_tvec)
        # self.camLight.setIntensity(100)
        # for mol in mols:
        self.setRootEntity(self.rootEntity)

        # Capture entity for screenshot
        # self.__capture = Qt3DRender.QRenderCapture(self.__cam)
        self.__capture = Qt3DRender.QRenderCapture(self.__cam)
        self.activeFrameGraph().setParent(self.__capture)
        self.setActiveFrameGraph(self.__capture)

        # Add shortcuts
        self.keyExport = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_E), self)
        self.keyExport.activated.connect(self.__export_pov)
        self.keyPrint = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_P), self)
        self.keyPrint.activated.connect(self.__capture_scene)

    def add_vibmode(self, vib_mode: Type1Vib,
                    model: tp.Optional[str] = None,
                    material: tp.Optional[str] = None,
                    color: tp.Optional[tp.Any] = None,
                    color2: tp.Optional[tp.Any] = None,
                    mol: int = -1):
        """Add vibrational mode.

        Creates a VibMode class and populate.

        Parmaters
        ---------
        vib_mode
            List of atomic displacements for a given mode.
        model
            Representation model to visualize the modes.
        material
            Material to use in the representation of the modes.
        color
            Color of the object representing the displacement.
        color2
            Secondary color of the object representing the displacement.
        mol
            Set target molecule if several of interest.
            Default: last included molecule.
        """
        if self.__vib is not None:
            del self.__vib
        if self.__nmols > 1:
            try:
                atcrd = self.__mol[mol].at_crd
            except IndexError as err:
                raise ValueError('Molecule does not exist.') from err
        else:
            atcrd = self.__mol.at_crd
        self.__vib = VibMode(atcrd, vib_mode, model, material,
                             color=color, color2=color2,
                             rootEntity=self.rootEntity)

    def upd_vibmode(self, vib_mode: Type1Vib):
        if self.__vib is None:
            self.__add_vibmode(vib_mode)
        else:
            self.__vib.update_vib(vib_mode)

    @QtCore.Slot(list)
    def atom_clicked(self, molatom: tp.Sequence[int]):
        """Capture the information on the clicked atom.

        Captures the information from the clicked atom on one of the
        molecules.

        Parameters
        ----------
        molatom
            Molar atom, as a list of integers.
        """
        self.click_mol.emit(molatom)

    def mousePressEvent(self, mouseEvent):
        """Define the event for mouse pressing."""
        if mouseEvent.button() == QtCore.Qt.RightButton:
            self.__cam.setViewCenter(QtGui.QVector3D(0, 0, 0))
        # print(mouseEvent.x())
        # self.__cam.setViewCenter(QtGui.QVector3D(mouseEvent.x(),
        #                                             mouseEvent.y(), 0))
        super(MolWin, self).mousePressEvent(mouseEvent)

    def __export_pov(self):
        """Export to POV-Ray description file.

        Builds a POV-Ray description file based on current
        configuration.
        """
        workdir = os.getcwd()
        if self.__nmols > 1:
            povfile = 'scene.pov'
        elif self.__fnames is None:
            povfile = 'molecule.pov'
        else:
            povfile = os.path.splitext(self.__fnames)[0] + '.pov'
            workdir = os.path.realpath(os.path.dirname(povfile))

        dialog = QtWidgets.QFileDialog()
        dialog.setDirectory(workdir)
        dialog.setWindowTitle('Export POV-Ray description file')
        dialog.setNameFilter("POV-Ray files (*.pov)")
        dialog.setNameFilter("All files (*)")
        dialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptMode.AcceptSave)
        dialog.setDefaultSuffix('.pov')
        dialog.selectFile(os.path.basename(povfile))
        if not dialog.exec_():
            return

        povfile = os.path.realpath(dialog.selectedFiles()[0])

        rmat = self.__cam.transform().rotation().toRotationMatrix().data()
        # We need to transpose the rotation matrix since it is intended for
        # the camera and we want to apply it to the system.
        # We then permute Y and Z in each dimension to correct from the right
        # to the left-handed system.
        rotmat = [
            [rmat[0], rmat[6], rmat[3]],
            [rmat[2], rmat[8], rmat[5]],
            [rmat[1], rmat[7], rmat[4]]
        ]
        pov_opts = {
            'rotation': rotmat,
            'verbose': False
        }
        if self.__vib is not None:
            mode = self.__vib.mode
            pov_opts['id_vib'] = 0
            pov_opts['vib_model'] = self.__vib.viz_options('model') or 'arrows'
            pov_opts['vib_mater'] = self.__vib.viz_options('material') \
                or 'plastic'
        else:
            mode = None
        if self.__nmols == 1:
            pov_opts['mol_model'] = self.__mol.viz_options('model')
            pov_opts['mol_mater'] = self.__mol.viz_options('material')
            builder = POVBuilder(None, povfile=povfile,
                                 at_crds=self.__mol.at_crd,
                                 at_labs=self.__mol.at_lab,
                                 v_modes=mode,
                                 nmols=self.__nmols, in_au=False)
            builder.write_pov(**pov_opts)
        else:
            pov_opts['mol_model'] = self.__mol[0].viz_options('model')
            pov_opts['mol_mater'] = self.__mol[0].viz_options('material')
            builder = POVBuilder(None, povfile=povfile,
                                 at_crds=[item.at_crd for item in self.__mol],
                                 at_labs=[item.at_lab for item in self.__mol],
                                 v_modes=mode,
                                 nmols=self.__nmols, in_au=False)
            builder.write_pov(**pov_opts)

        txt = MSG_POV_CMD_HTML.format(file=povfile)
        if os.sys.platform == 'win32':
            txt += '\n\n' + MSG_POV_WIN
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Information)
        msgBox.setText("Information for the execution of POV-Ray")
        msgBox.setInformativeText(txt)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.Ok)
        msgBox.exec()

    def __capture_scene(self):
        """Capture current scene."""
        shot = self.__capture.requestCapture()
        loop = QtCore.QEventLoop()
        shot.completed.connect(loop.quit)
        QtCore.QTimer.singleShot(200, loop.quit)  # Timeout after 1 second
        loop.exec()

        img = shot.image()
        if not img:
            return

        supported = ('.bmp', '.jpeg', '.jpg', '.png', '.ppm', '.xbm', '.xpm')
        ftype = 'png'
        if self.__nmols > 1 or self.__fnames is None:
            fname = f'scene.{ftype}'
        else:
            fname = f'{os.path.splitext(self.__fnames)[0]}.{ftype}'

        dialog = QtWidgets.QFileDialog()
        dialog.setDirectory(os.getcwd())
        dialog.setWindowTitle('Export screenshot')
        dialog.setNameFilter(
            f"Supported image files (*{' *'.join(supported)})")
        dialog.setNameFilter("All files (*)")
        dialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptMode.AcceptSave)
        dialog.setDefaultSuffix('.png')
        dialog.selectFile(os.path.basename(fname))
        if not dialog.exec_():
            return

        fname = os.path.realpath(dialog.selectedFiles()[0])
        img.save(fname)
