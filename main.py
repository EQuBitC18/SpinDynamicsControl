# hello_world.py

from PyQt6.QtWidgets import QApplication, QWidget


# create the QApplication
app = QApplication([])

# create the main window
window = QWidget(windowTitle='Hello World')
window.show()

# start the event loop
app.exec()















