# from backend import *
# import os
from ui import *
from PyQt5.QtWidgets import QApplication

def main():
    app = QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()