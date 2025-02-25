"""Build user interfaces for molecular visualization.

Provides classes to build simple or more advanced interfaces to
visualize molecules and associated properties.
"""
import os
import typing as tp

import numpy as np

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DRender import Qt3DRender
from PySide6.Qt3DExtras import Qt3DExtras

from estampes.base.types import Type1Vib, TypeAtCrdM, TypeAtLabM, TypeBondsM, \
    TypeColor
from estampes.tools.atom import convert_labsymb
from estampes.tools.math import rotate
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
            `skip_guide`: Skip introduction message.
        """
        super(MolWin, self).__init__()

        self.__vib = None
        self.__mol = None
        self.__nmols = 0
        self.__vib_changed_since_anim = None
        self.__dialog = None

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
        self.__cam_origin = np.array([0, 0, 40])
        self.__cam = self.camera()
        self.__cam.lens().setPerspectiveProjection(45, 16 / 9, 0.1, 1000)
        self.__cam.setPosition(QtGui.QVector3D(*self.__cam_origin))
        self.__cam.setUpVector(QtGui.QVector3D(0, 1, 0))
        self.__cam.setViewCenter(QtGui.QVector3D(0, 0, 0))

        # For camera controls
        self.rootEntity = Qt3DCore.QEntity()
        Molecule.reset_counter()
        if nmols == 1:
            self.__nmols = 1
            self.__mol = Molecule(atlabs, atcrds, bonds, model, material,
                                  col_bond_as_atom, molcols, self.rootEntity)
            self.__mol.addMouse(self.__cam)
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
        self.light_transfo = Qt3DCore.QTransform()
        self.light_transfo.setTranslation(QtGui.QVector3D(0, 50, 100))
        self.obj_light.addComponent(self.camLight)
        self.obj_light.addComponent(self.light_transfo)
        # self.camLight.setIntensity(100)
        # for mol in mols:
        self.setRootEntity(self.rootEntity)

        # Capture entity for screenshot
        # self.__capture = Qt3DRender.QRenderCapture(self.__cam)
        self.__capture = Qt3DRender.QRenderCapture(self.__cam)
        self.activeFrameGraph().setParent(self.__capture)
        self.setActiveFrameGraph(self.__capture)

        if not extras.get('skip_guide', False):
            self.__show_help()

        # Add shortcuts
        self.keyHelp = QtGui.QShortcut(
            QtGui.QKeySequence(
                QtGui.Qt.CTRL | QtGui.Qt.SHIFT | QtGui.Qt.Key_H),
            self,
            context=QtCore.Qt.ApplicationShortcut)
        self.keyHelp.activated.connect(self.__show_help)
        self.keyHelp.activatedAmbiguously.connect(self.__show_help)
        self.keyVib = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_A), self,
            context=QtCore.Qt.ApplicationShortcut)
        self.keyVib.activated.connect(self.__anim_vibration)
        self.keyExport = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_E), self,
            context=QtCore.Qt.ApplicationShortcut)
        self.keyExport.activated.connect(self.__export_pov)
        self.keyGeom = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_G), self,
            context=QtCore.Qt.ApplicationShortcut)
        self.keyGeom.activated.connect(self.__show_geom)
        self.keyPrint = QtGui.QShortcut(
            QtGui.QKeySequence(QtGui.Qt.CTRL | QtGui.Qt.Key_P), self,
            context=QtCore.Qt.ApplicationShortcut)
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
        self.__vib_changed_since_anim = True
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

    def upd_vibmode(self, vib_mode: Type1Vib, **kwargs):
        """Update vibrational mode.

        Extra keywords are sent to add_vibmode.
        """
        self.__vib_changed_since_anim = True
        if self.__vib is None:
            self.add_vibmode(vib_mode, **kwargs)
        else:
            self.__vib.update_vib(vib_mode)
            if self.__mol.animation_status != 'unknown':
                self.__mol.redraw()
        if self.__mol.animation_status == 'running':
            self.__mol.animation_status = 'reset'
            self.__anim_vibration()
        elif self.__mol.animation_status != 'unknown':
            self.__mol.animation_status = 'reset'

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

        self.__dialog = QtWidgets.QFileDialog()
        self.__dialog.setDirectory(workdir)
        self.__dialog.setWindowTitle('Export POV-Ray description file')
        self.__dialog.setNameFilter("POV-Ray files (*.pov)")
        self.__dialog.setNameFilter("All files (*)")
        self.__dialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        self.__dialog.setAcceptMode(
            QtWidgets.QFileDialog.AcceptMode.AcceptSave)
        self.__dialog.setDefaultSuffix('.pov')
        self.__dialog.selectFile(os.path.basename(povfile))
        if not self.__dialog.exec_():
            return

        povfile = os.path.realpath(self.__dialog.selectedFiles()[0])

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

        self.__dialog = QtWidgets.QFileDialog()
        self.__dialog.setDirectory(os.getcwd())
        self.__dialog.setWindowTitle('Export screenshot')
        self.__dialog.setNameFilter(
            f"Supported image files (*{' *'.join(supported)})")
        self.__dialog.setNameFilter("All files (*)")
        self.__dialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        self.__dialog.setAcceptMode(
            QtWidgets.QFileDialog.AcceptMode.AcceptSave)
        self.__dialog.setDefaultSuffix('.png')
        self.__dialog.selectFile(os.path.basename(fname))
        if not self.__dialog.exec():
            return

        fname = os.path.realpath(self.__dialog.selectedFiles()[0])
        img.save(fname)

    def __anim_vibration(self):
        """Animate vibration.

        Runs/stop the animation of the vibration
        """
        if self.__vib is None:
            return
        if self.__nmols == 1:
            if self.__vib_changed_since_anim:
                self.__mol.animate(self.__vib.mode)
                self.__vib_changed_since_anim = False
            else:
                self.__mol.animate()
        else:
            raise NotImplementedError('Multi-mols vibrations NYI')

    def __show_geom(self):
        """Open a dialog displaying the current geometry."""
        def update_geom(which: int = 0):
            if 0 <= which < 3:
                atcrd = self.__mol.get_at_crd(True)
            else:
                atcrd = self.__mol.get_at_crd(False)
            if which % 3 != 2:
                transformation = self.__cam.transform()
                rmat = np.reshape(
                    transformation.rotation().toRotationMatrix().data(),
                    (3, 3))
                atcrd = rotate(rmat, atcrd)
                if which % 3 == 0:
                    tvec = np.array([
                        0, 0,
                        transformation.translation().z()-self.__cam_origin[-1]
                        ])
                    atcrd += tvec
            lines = [f'<pre>{len(atlab)}\n', f'{list_choices[which]}\n']
            for lab, crd in zip(atlab, atcrd):
                lines.append(fmt_xyz.format(at=lab, xyz=crd))
            lines.append('</pre>')
            geom.setText(''.join(lines))

        def to_clipboard():
            clipboard = QtGui.QClipboard()
            clipboard.setText(geom.toPlainText())

        def save_xyz():
            if self.__nmols > 1 or self.__fnames is None:
                fname = 'system.xyz'
            else:
                fname = f'{os.path.splitext(self.__fnames)[0]}.xyz'
            win_save = QtWidgets.QFileDialog()
            win_save.setDirectory(os.getcwd())
            win_save.setWindowTitle("Export XYZ file")
            win_save.setNameFilter("XYZ files (*.xyz *.XYZ)")
            win_save.setNameFilter("All files (*)")
            win_save.setFileMode(QtWidgets.QFileDialog.AnyFile)
            win_save.setAcceptMode(QtWidgets.QFileDialog.AcceptMode.AcceptSave)
            win_save.setDefaultSuffix('.xyz')
            win_save.selectFile(os.path.basename(fname))
            if not win_save.exec_():
                return
            fname = os.path.realpath(win_save.selectedFiles()[0])
            with open(fname, 'w', encoding='utf-8') as fobj:
                fobj.write(geom.toPlainText())

        if self.__nmols > 1:
            raise NotImplementedError('Geometries for multiple molecules NYI')

        fmt_xyz = '{at:<2s}{xyz[0]:14.8f}{xyz[1]:14.8f}{xyz[2]:14.8f}\n'
        list_choices = [
            'Current geometry, including all transformations',
            'Current geometry, including only rotations',
            'Current geometry in original orientation',
            'Reference geometry, including all transformations',
            'Reference geometry, including only rotations',
            'Reference geometry in original orientation'
        ]
        atlab = self.__mol.at_lab
        for atom in atlab:
            if isinstance(atom, int):
                do_convert = True
                break
        else:
            do_convert = False
        if do_convert:
            atlab = convert_labsymb(True, *atlab)
        self.__dialog = QtWidgets.QDialog()
        layout = QtWidgets.QVBoxLayout()
        geom = QtWidgets.QTextEdit()
        geom.setReadOnly(True)
        update_geom(0)
        layout.addWidget(geom)

        choices = QtWidgets.QComboBox()
        choices.addItems(list_choices)
        choices.currentIndexChanged.connect(update_geom)
        layout.addWidget(choices)

        buttons = QtWidgets.QHBoxLayout()
        btn_close = QtWidgets.QPushButton('Close')
        btn_close.clicked.connect(self.__dialog.close)
        btn_clip = QtWidgets.QPushButton('Copy to clipboard')
        btn_clip.clicked.connect(to_clipboard)
        btn_save = QtWidgets.QPushButton('Save')
        btn_save.clicked.connect(save_xyz)
        buttons.addWidget(btn_close)
        buttons.addStretch()
        buttons.addWidget(btn_clip)
        buttons.addStretch()
        buttons.addWidget(btn_save)
        layout.addLayout(buttons)

        self.__dialog.setLayout(layout)
        self.__dialog.exec()

    def __show_help(self):
        """Show help message"""
        self.__dialog = QtWidgets.QMessageBox()
        self.__dialog.setIcon(QtWidgets.QMessageBox.Information)
        self.__dialog.setText("Available commands")
        self.__dialog.setInformativeText("""
<style>
h3 {
    font-style: italic;
    text-decoration: underline double;
}
dt {
    font-weight: bold;
    font-size: 1.2em;
    margin-bottom: .2em;
    margin-top: .2em;
}
</style>
<h3>Basic navigation</h3>

<dl>
    <dt>Left mouse button</dt>
    <dd>While the left mouse button is pressed, mouse movement along x-axis \
moves the camera left and right and movement along y-axis moves it up and \
down.</dd>
    <dt>Right mouse button</dt>
    <dd>While the right mouse button is pressed, mouse movement along x-axis \
pans the camera around the camera view center and movement along y-axis tilts \
it around the camera view center.</dd>
    <dt>Both left and right mouse button</dt>
    <dd>While both the left and the right mouse button are pressed, mouse \
movement along y-axis zooms the camera in and out without changing the view \
center.</dd>
    <dt>Mouse scroll wheel</dt>
    <dd>Zooms the camera in and out without changing the view center.</dd>
    <dt>Arrow keys</dt>
    <dd>Move the camera vertically and horizontally relative to camera \
viewport.</dd>
    <dt>Page up and page down keys</dt>
    <dd>Move the camera forwards and backwards.</dd>
    <dt>Shift key</dt>
    <dd>Changes the behavior of the up and down arrow keys to zoom the camera \
in and out without changing the view center. The other movement keys are \
disabled.</dd>
    <dt>Alt key</dt>
    <dd>Changes the behovior of the arrow keys to pan and tilt the camera \
around the view center. Disables the page up and page down keys.</dd>
    <dt>Escape</dt>
    <dd>Moves the camera so that entire scene is visible in the camera \
viewport.</dd>
</dl>

<h3>Shortcuts</h3>

<p>On Mac OSX, <code>Ctrl</code> is replaced by <code>Command</code>.</p>

<dl>
    <dt>Ctrl+A</dt>
    <dd>Animates vibration if a vibration has been selected.<br />
Pressing it again pauses or resumes the motion.</dd>
    <dt>Ctrl+E</dt>
    <dd>Exports current scene as a POV-Ray scene description file.</dd>
    <dt>Ctrl+G</dt>
    <dd>Displays the current geometry as a text file.</dd>
    <dt>Ctrl+Shift+H</dt>
    <dd>Prints this help message.</dd>
    <dt>Ctrl+P</dt>
    <dd>Captures the current scene and saves it as an image.<br />
    <strong>Note:</strong> The option is currently experimental and may \
fail.</dd>
</dl>
""")
        self.__dialog.setStandardButtons(QtWidgets.QMessageBox.Close)
        self.__dialog.setDefaultButton(QtWidgets.QMessageBox.Close)
        self.__dialog.show()
