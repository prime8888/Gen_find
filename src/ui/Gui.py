import sys
import csv
from ui.utils import toggle_all_checkboxes, checkbox_toggled
# import requests
from PyQt5.QtWidgets import QApplication, QMainWindow, QTreeWidget, QPushButton, QMenu, QTreeWidgetItem, QVBoxLayout, QWidget, QTextEdit, QHBoxLayout, QLabel, QCheckBox, QFrame, QSizePolicy, QMessageBox, QTreeView, QFileSystemModel
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QDir
from backend import *
# from Bio import Entrez

# Define your email here for NCBI access
# Entrez.email = "emailbidon@ncbi.nlm.nih.gov"

# class FetchThread(QThread):
#     dataFetched = pyqtSignal(dict)
#     errorOccurred = pyqtSignal(str)

#     def __init__(self, tax_ids):
#         super().__init__()
#         self.tax_ids = tax_ids

#     def run(self):
#         try:
#             result = {}
#             for tax_name, tax_id in self.tax_ids.items():
#                 details = self.fetch_tax_details(tax_id)
#                 result[tax_name] = details
#             self.dataFetched.emit(result)
#         except Exception as e:
#             self.errorOccurred.emit(str(e))

#     def fetch_tax_details(self, tax_id, details=None):
#         if details is None:
#             details = []
#         handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
#         records = Entrez.read(handle)
#         handle.close()
#         for record in records:
#             if 'LineageEx' in record:
#                 for lineage in record['LineageEx']:
#                     details.append((lineage['ScientificName'], lineage['TaxId'], lineage['Rank']))
#                     if lineage['Rank'] in ['family', 'genus']:  # You can adjust the depth condition
#                         self.fetch_tax_details(lineage['TaxId'], details)
#         return details

class WorkerThread(QThread):
    dataFetched = pyqtSignal()
    errorOccurred = pyqtSignal(str)

    def __init__(self, max_workers, selected_paths, selected_regions, rettype="gbwithparts"):
        super().__init__()
        self.max_workers = max_workers
        self.selected_paths = selected_paths
        self.rettype = rettype
        self.selected_regions = selected_regions

    def run(self):
        try:
            organisms = get_organisms_from_path_list(self.selected_paths)
            print(organisms)
            print(self.selected_regions)
            organisms = fetch_nc_ids_for_organisms(organisms)
            records = fetch_genbank_records(organisms, rettype=self.rettype, max_workers=self.max_workers)
            process_genomic_data(records, self.selected_regions)
            self.dataFetched.emit()
        except Exception as e:
            self.errorOccurred.emit(str(e))

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
        self.max_workers = os.cpu_count() + 2
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
        # Set up the file system model
        self.model = QFileSystemModel()
        self.model.setRootPath(QDir.rootPath())

        tree_layout = QVBoxLayout()
        # Set up the tree view
        self.tree = QTreeView()
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(os.path.join(QDir.currentPath(), "Results")))
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
        start_button = QPushButton("Start Fetch")
        stop_button = QPushButton("Stop Fetch")
        expand_button = QPushButton("Expand All")
        collapse_button = QPushButton("Collapse All")

        start_button.clicked.connect(self.start_clicked)
        stop_button.clicked.connect(self.stop_clicked)
        expand_button.clicked.connect(lambda: self.tree.expandAll())
        collapse_button.clicked.connect(lambda: self.tree.collapseAll())

        button_style = "QPushButton { background-color: #313042; color: white; font-weight: bold; }"
        for button in [start_button, stop_button, expand_button, collapse_button]:
            button.setStyleSheet(button_style)
            button_layout.addWidget(button)

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
        # item = self.tree.currentItem()
        # if item is not None:
        #     self.info_output.setText(f"You clicked on {item.text(0)}")
        indexes = self.tree.selectionModel().selectedIndexes()
        # Filter indexes to get only those from the first column
        unique_indexes = {idx.row(): idx for idx in indexes if idx.column() == 0}
        paths = [self.model.filePath(idx) for idx in unique_indexes.values()]
        self.selected_paths = paths
        print(paths)

    def start_clicked(self):
        if not self.selected_paths:
            QMessageBox.critical(self, "Error", "Please select a path")
            return
        if not any(checkbox.isChecked() for checkbox in self.checkboxes):
            QMessageBox.critical(self, "Error", "Please select at least one region")
            return
        self.selected_regions = [checkbox.text() for checkbox in self.checkboxes if checkbox.isChecked()]
        if "ALL" in self.selected_regions:
            self.selected_regions = ["mobile_element", "5'UTR", "telomere", "intron", "3'UTR", "rRNA", "centromere", "tRNA", "ncRNA", "CDS"]
        self.thread = WorkerThread(self.max_workers, self.selected_paths, self.selected_regions)
        self.thread.errorOccurred.connect(self.handle_error)
        self.thread.start()

    def stop_clicked(self):
        if self.thread.isRunning():
            self.thread.terminate()
        print("Fetch stopped")

    def update_tree(self, fetched_data):
        self.tree.clear()
        for tax_name, details in fetched_data.items():
            parent_item = QTreeWidgetItem(self.tree)
            parent_item.setText(0, tax_name)
            self.add_tax_items(parent_item, details)

    def add_tax_items(self, parent_item, details):
        for name, tax_id, rank in details:
            item = QTreeWidgetItem(parent_item)
            item.setText(0, f"{name} ({rank})")

    def add_tree_item(self, parent_item, text):
        item = QTreeWidgetItem(parent_item)
        item.setText(0, text)

    def handle_error(self, error_message):
        QMessageBox.critical(self, "Error", f"An error occurred: {error_message}")
        self.log_text.append(f"Error occurred: {error_message}")

    def create_information_section(self, layout):
        info_frame = QFrame()
        info_frame.setStyleSheet("border: none;")
        info_layout = QVBoxLayout(info_frame)

        info_title = QLabel("Information")
        info_title.setObjectName("info_title")
        info_title.setAlignment(Qt.AlignCenter)
        info_title.setStyleSheet("padding: 6px; border: none;")
        self.info_output = QTextEdit()
        self.info_output.setStyleSheet("background-color:#1a1a2e;")
        self.info_output.setReadOnly(True)

        info_layout.addWidget(info_title)
        info_layout.addWidget(self.info_output)
        layout.addWidget(info_frame, 1)

    def create_logs_section(self, layout):
        log_frame = QFrame()
        log_frame.setStyleSheet("border: none;")
        log_vbox = QVBoxLayout(log_frame)
    
        log_title = QLabel("Logs")
        log_title.setObjectName("log_title")
        log_title.setAlignment(Qt.AlignCenter)
        log_title.setStyleSheet("padding: 6px;")
        self.log_text = QTextEdit()  # Change here to use self.log_text
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("background-color:#1a1a2e;")
        log_vbox.addWidget(log_title)
        log_vbox.addWidget(self.log_text)
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
                checkbox.stateChanged.connect(lambda state, chks=self.checkboxes[:-1], log_text=self.log_text: toggle_all_checkboxes(chks, state, log_text))
            else:
                self.checkboxes.append(checkbox)
                checkbox.stateChanged.connect(lambda state, opt=option, chk=checkbox, all_chk=self.checkboxes[-1], log_text=self.log_text:
                                              checkbox_toggled(opt, state, chk, all_chk, log_text))
            checkbox_layout.addWidget(checkbox)
    
        selection_vbox.addWidget(checkbox_frame)
        checkbox_frame.setMinimumHeight(400)
        checkbox_frame.setStyleSheet(" background-color:#1a1a2e ;") 
        checkbox_frame.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        layout.addWidget(selection_frame, 1)

    def apply_stylesheet(self):
        stylesheet_path = "./ui/style.qss"
        with open(stylesheet_path, "r") as file:
            self.setStyleSheet(file.read())

