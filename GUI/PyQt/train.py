# Simple hand made app without designer

from PyQt5 import QtWidgets  # to add widgets
from PyQt5.QtWidgets import QApplication, QMainWindow  # to create App itself and window
#  pyuic5 name.ui -o name.py - запускаем из папки с файлом ui в conda prompt (if anaconda is not installed then use cmd)

import sys

class Window(QMainWindow):
    def __init__(self):
        super(Window, self).__init__()

        self.setWindowTitle('Best name')
        self.setGeometry(300, 250, 300, 250)  # (pixel offset from top left courner to the right, to the bottom, size width, size height)

        self.new_text = QtWidgets.QLabel(self)

        # widgets must be added after the 'window' object and before window.show() method
        self.main_text = QtWidgets.QLabel(self)
        self.main_text.setText('label text')
        self.main_text.move(100, 100)
        self.main_text.adjustSize()

        self.btn = QtWidgets.QPushButton(self)
        self.btn.move(100, 0)
        self.btn.setText('Push')
        self.btn.clicked.connect(self.button1)

    def button1(self):
        self.new_text.setText('HERE')
        self.new_text.move(10, 10)
        self.new_text.adjustSize()

def application():
    app = QApplication(sys.argv)  # create object of our app and pass system data to it (info about computer)
    window = Window()

    window.show()
    sys.exit(app.exec_())


# if __name__ == "__main__":
#     application()

import pyqtgraph.examples
pyqtgraph.examples.run()

