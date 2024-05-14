import sys
import csv
from ui.utils import toggle_all_checkboxes, checkbox_toggled
# import requests
from PyQt5.QtWidgets import QMainWindow, QPushButton, QMenu, QTreeWidgetItem, QVBoxLayout, QWidget, QTextEdit, QHBoxLayout, QLabel, QCheckBox, QFrame, QSizePolicy, QMessageBox, QTreeView, QFileSystemModel, QProgressBar
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QDir, QRunnable, QThreadPool, QMetaObject, QObject, Q_ARG
from backend import *
from ui.logger import *


class Signals(QObject):
    processComplete = pyqtSignal(object, object)  # You can customize these signals as needed
    processFailed = pyqtSignal(str)

class ProcessRunnable(QRunnable):
    def __init__(self, organism, id, path, selected_regions, rettype):
        super(ProcessRunnable, self).__init__()
        self.signals = Signals()
        self.organism = organism
        self.id = id
        self.path = path
        self.rettype = rettype
        self.selected_regions = selected_regions
        self._is_running = True

    def run(self):
        if not self._is_running:
            return
        
        try:
            fetch_records_and_process(self.organism, self.id, self.path, self.selected_regions, rettype=self.rettype)
            self.signals.processComplete.emit(self.organism, self.path)
        except Exception as e:
            self.signals.processFailed.emit(str(e))

    def stop(self):
        self._is_running = False


class WorkerThread(QThread):
    dataFetched = pyqtSignal()
    errorOccurred = pyqtSignal(str)
    fullStop = pyqtSignal()
    progressUpdated = pyqtSignal(int)

    def __init__(self, selected_paths, selected_regions, rettype="gbwithparts", max_workers=5):
        super(WorkerThread, self).__init__()
        self.selected_paths = selected_paths
        self.selected_regions = selected_regions
        self.rettype = rettype
        self.threadPool = QThreadPool.globalInstance()
        self.threadPool.setMaxThreadCount(max_workers)
        self.runnables = []
        self.logger = logging.getLogger()
        self.stopped = False
        self.totalTasks = 0
        self.completedTasks = 0

    def run(self):
        try:
            organisms = get_organisms_from_path_list(self.selected_paths)
            organisms = fetch_nc_ids_for_organisms(organisms)
            self.totalTasks = sum(len(data['ids']) for data in organisms.values())

            for organism, data in organisms.items():
                for id in data['ids']:
                    runnable = ProcessRunnable(organism, id, data['path'], self.selected_regions, self.rettype)
                    runnable.signals.processComplete.connect(self.recordProcessed)
                    runnable.signals.processFailed.connect(self.handleError)
                    self.runnables.append(runnable)
                    self.threadPool.start(runnable)
            self.threadPool.waitForDone()
            if not self.stopped:
                self.dataFetched.emit()
            else:
                self.fullStop.emit()
            
        except Exception as e:
            self.errorOccurred.emit(str(e))

    def recordProcessed(self, organism, path):
        self.completedTasks += 1
        progress = int((self.completedTasks / self.totalTasks) * 100)
        self.progressUpdated.emit(progress)

    def handleError(self, error):
        self.errorOccurred.emit(error)

    def stop(self):
        self.stopped = True
        self.logger.warning("Stopping all tasks. Waiting for working threads to finish...")
        for runnable in self.runnables:
            if hasattr(runnable, 'stop'):
                runnable.stop()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # self.tax_ids = {'Humans': '9606', 'Fruit flies': '7227'}  # Example tax IDs
        self.initUI()
        self.onStartApp()

    def initUI(self):
        self.setup_central_widget()
        self.setGeometry(300, 200, 1100, 800)
        self.setWindowTitle('GENOME')
        self.apply_stylesheet()
        self.show()

    def onStartApp(self):
        self.selected_regions = []
        self.selected_paths = []
        # self.max_workers = os.cpu_count() + 8
        self.max_workers = 12
        self.logger = logging.getLogger()
        fetch_and_update_overview()


    def setup_central_widget(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        tree_layout = self.setup_tree_layout()
        main_layout.addLayout(tree_layout, 1)

        right_side_layout = QVBoxLayout()
        self.create_information_section(right_side_layout)
        logs_selection_layout = QHBoxLayout()
        self.create_logs_section(logs_selection_layout)
        self.create_selection_section(logs_selection_layout)
        right_side_layout.addLayout(logs_selection_layout)
        main_layout.addLayout(right_side_layout, 2)

    def setup_tree_layout(self):
        start_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        if start_dir.endswith('src'):
            # Move up to the parent directory
            start_dir = os.path.dirname(start_dir)
        else:
            start_dir = start_dir
        # Set up the file system model
        start_dir = os.path.join(start_dir, "Results")
        self.model = QFileSystemModel()
        self.model.setRootPath(start_dir)

        tree_layout = QVBoxLayout()
        # Set up the tree view
        self.tree = QTreeView()
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(start_dir))
        for column in range(1, self.model.columnCount()):
            self.tree.hideColumn(column)
        # self.tree.setHeaderLabel("Taxonomy Hierarchy")
        self.tree.setSelectionMode(QTreeView.ExtendedSelection)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self.open_menu)
        self.tree.clicked.connect(self.on_tree_clicked)
        self.tree.setHeaderHidden(False)
        self.tree.setStyleSheet("background-color:#1a1a2e;")
        self.tree.header().setStyleSheet("QHeaderView::section { background-color: #313042; border: none; font-size: 18px; font-weight: bold; }")
        self.tree.setRootIsDecorated(True)
        self.tree.setIndentation(20)
        tree_layout.addWidget(self.tree)

        button_layout = self.setup_buttons()
        tree_layout.addLayout(button_layout)
        return tree_layout

    def setup_buttons(self):
        button_layout = QHBoxLayout()
        start_button = QPushButton("Start Processing")
        self.start_button = start_button

        start_button.clicked.connect(self.start_clicked)

        button_style = "QPushButton { background-color: #313042; color: white; font-weight: bold; }"
        start_button.setStyleSheet(button_style)
        button_layout.addWidget(start_button)

        return button_layout

    def open_menu(self, position):
        menu = QMenu()
        expand_action = menu.addAction("Expand")
        collapse_action = menu.addAction("Collapse")
        action = menu.exec_(self.tree.viewport().mapToGlobal(position))
        if action == expand_action:
            self.tree.expandItem(self.tree.currentItem())
        elif action == collapse_action:
            self.tree.collapseItem(self.tree.currentItem())

    def on_tree_clicked(self, index):
        indexes = self.tree.selectionModel().selectedIndexes()
        # Filter indexes to get only those from the first column
        unique_indexes = {idx.row(): idx for idx in indexes if idx.column() == 0}
        paths = [self.model.filePath(idx) for idx in unique_indexes.values()]
        self.selected_paths = paths
        self.logger.info(f"Selected paths:")
        for path in paths:
            self.logger.info(f"{path}")

    def start_clicked(self):
        if not self.selected_paths:
            QMessageBox.critical(self, "Error", "Please select a path")
            return
        if not any(checkbox.isChecked() for checkbox in self.checkboxes):
            QMessageBox.critical(self, "Error", "Please select at least one region")
            return
        self.start_button.setText("Stop Processing")
        try:        
            self.start_button.clicked.disconnect()
        except TypeError:
            pass
        self.start_button.clicked.connect(self.stop_clicked)
        self.selected_regions = [checkbox.text() for checkbox in self.checkboxes if checkbox.isChecked()]
        if "ALL" in self.selected_regions:
            self.selected_regions = ["mobile_element", "5'UTR", "telomere", "intron", "3'UTR", "rRNA", "centromere", "tRNA", "ncRNA", "CDS"]
        self.thread = WorkerThread(self.selected_paths, self.selected_regions, max_workers=self.max_workers)
        self.thread.errorOccurred.connect(self.handle_error)
        self.thread.dataFetched.connect(self.processing_finished)
        self.thread.fullStop.connect(self.processing_stopped)
        self.thread.progressUpdated.connect(self.progressBar.setValue)
        self.thread.start()

    def processing_finished(self):
        self.progressBar.setValue(100)
        self.logger.info("All records parsed successfully.")
        self.start_button.setText("Start Processing")
        try:        
            self.start_button.clicked.disconnect()
        except TypeError:
            pass
        self.start_button.clicked.connect(self.start_clicked)
        self.start_button.setDisabled(False)

    def processing_stopped(self):
        self.logger.info("All remaining threads stopped.")
        self.start_button.setText("Start Processing")

        try:        
            self.start_button.clicked.disconnect()
        except TypeError:
            pass
        
        self.start_button.clicked.connect(self.start_clicked)
        self.start_button.setDisabled(False)

    def stop_clicked(self):
        if self.thread:
            self.thread.stop()
            self.start_button.setText("Waiting for remaining threads to stop...")
            self.start_button.setDisabled(True)
            try:        
                self.start_button.clicked.disconnect()
            except TypeError:
                pass



    def handle_error(self, error_message):
        self.progressBar.reset()
        QMessageBox.critical(self, "Error", f"An error occurred: {error_message}")
        self.log_text.append(f"Error occurred: {error_message}")
        self.start_button.setText("Start Processing")
        try:        
            self.start_button.clicked.disconnect()
        except TypeError:
            pass
        self.start_button.clicked.connect(self.start_clicked)
        self.start_button.setDisabled(False)

    def create_information_section(self, layout):
        info_frame = QFrame()
        info_frame.setStyleSheet("border: none;")
        info_layout = QVBoxLayout(info_frame)

        # Title Label
        info_title = QLabel("Information")
        info_title.setObjectName("info_title")
        info_title.setAlignment(Qt.AlignCenter)
        info_title.setStyleSheet("padding: 6px; border: none; font-weight: bold;")

        # Progress Bar setup
        self.progressBar = QProgressBar()
        self.progressBar.setMaximumWidth(250)  # Limit the width
        self.progressBar.setMinimumHeight(20)  # Ensure there's a minimum height
        self.progressBar.setMinimumSize(250, 20)  # Minimum size to ensure visibility
        self.progressBar.setMaximumSize(500, 30)  # Optional: limit the maximum size if necessary
        self.progressBar.setMaximum(100)
        self.progressBar.setValue(0)  # Set an initial value to make it visible for testing
        self.progressBar.setStyleSheet("""
            QProgressBar {
                border: 2px solid grey;
                border-radius: 5px;
                background-color: #f0f0f0; 
                color: black; 
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: green; /* Color of the progress meter */
                width: 20px; /* Chunk size */
            }
        """)

        # Horizontal layout to center the progress bar
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.progressBar)
        hbox.addStretch(1)

        # Adding widgets to the vertical layout of the frame
        info_layout.addWidget(info_title)
        info_layout.addLayout(hbox)  # Add the centered progress bar
        info_layout.addStretch(1)

        # Add the frame to the main layout
        layout.addWidget(info_frame)

    def create_logs_section(self, layout):
        log_frame = QFrame()
        log_frame.setStyleSheet("border: none;")
        log_vbox = QVBoxLayout(log_frame)
    
        log_title = QLabel("Logs")
        log_title.setObjectName("log_title")
        log_title.setAlignment(Qt.AlignCenter)
        log_title.setStyleSheet("padding: 6px;")
        self.log_text_edit = QTextEdit()
        self.log_text_edit.setReadOnly(True)
        logTextBox = QTextEditLogger(self.log_text_edit)
        
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        logTextBox.setFormatter(formatter)
        self.log_text_edit.setStyleSheet("background-color:#1a1a2e;")
        
        logger = logging.getLogger()
        logger.addHandler(logTextBox)
        logger.setLevel(logging.DEBUG)

        log_vbox.addWidget(log_title)
        log_vbox.addWidget(self.log_text_edit)
        layout.addWidget(log_frame, 1)
    
    def create_selection_section(self, layout):
        selection_frame = QFrame()
        selection_vbox = QVBoxLayout(selection_frame)
        selection_title = QLabel("Regions Fonctionnelles")
        selection_title.setObjectName("selection_title")
        selection_title.setAlignment(Qt.AlignCenter)
        selection_title.setStyleSheet("padding: 6px; border: none;")
        selection_vbox.addWidget(selection_title)

        checkbox_frame = QFrame()
        checkbox_layout = QVBoxLayout(checkbox_frame)
        options = ["mobile_element", "5'UTR", "telomere", "intron", "3'UTR", "rRNA", "centromere", "tRNA", "ncRNA", "CDS", "ALL"]
        self.checkboxes = []
    
        for option in options:
            checkbox = QCheckBox(option)
            checkbox.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
            if option == "ALL":
                checkbox.stateChanged.connect(lambda state, chks=self.checkboxes[:-1]: toggle_all_checkboxes(chks, state))
            else:
                self.checkboxes.append(checkbox)
                checkbox.stateChanged.connect(lambda state, opt=option, chk=checkbox, all_chk=self.checkboxes[-1]:
                                              checkbox_toggled(opt, state, chk, all_chk))
            checkbox_layout.addWidget(checkbox)
    
        selection_vbox.addWidget(checkbox_frame)
        checkbox_frame.setMinimumHeight(400)
        checkbox_frame.setStyleSheet(" background-color:#1a1a2e ;") 
        checkbox_frame.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        layout.addWidget(selection_frame, 1)

    def apply_stylesheet(self):
        # Load the stylesheet from the file
        project_path = os.path.dirname(os.path.dirname(__file__))
        stylesheet_path = os.path.join(project_path, "ui", "style.qss")
        with open(stylesheet_path, "r") as file:
            self.setStyleSheet(file.read())

