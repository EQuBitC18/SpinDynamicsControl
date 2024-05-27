import sys
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QGroupBox, QGridLayout, QComboBox, QLabel
)
from PyQt6.QtGui import QPixmap
from qiskit import QuantumCircuit
from qiskit.visualization import circuit_drawer


class QuantumCircuitGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Quantum Circuit GUI')
        
        # Set fixed size for the main window
        self.setFixedSize(1000, 800)
        
        # Main layout
        main_layout = QVBoxLayout()
        
        # Circuit Diagram Group
        self.circuit_diagram_group = QGroupBox('Circuit Diagram')
        self.circuit_diagram_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 16px;
                font-weight: bold;
            }
        """)
        circuit_diagram_layout = QVBoxLayout()
        
        self.circuit_diagram_label = QLabel('', self)
        circuit_diagram_layout.addWidget(self.circuit_diagram_label)
        
        self.circuit_diagram_group.setLayout(circuit_diagram_layout)
        main_layout.addWidget(self.circuit_diagram_group)
        
        # Bottom layout
        bottom_layout = QHBoxLayout()
        
        # Left Section of Bottom Layout
        left_bottom_layout = QVBoxLayout()
        
        # Parameter Configuration (Parent Group Box)
        param_config_group = QGroupBox('Parameter Configuration')
        param_config_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 16px;
                font-weight: bold;
            }
        """)
        param_config_layout = QHBoxLayout()

        # Hamiltonian Configuration (Subsection)
        hamiltonian_group = QGroupBox('Hamiltonian Configuration')
        hamiltonian_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 14px;
                font-weight: bold;
            }
        """)
        hamiltonian_layout = QGridLayout()
        hamiltonian_params = ['Interaction strengths:', 'Magnetic moments:', 'External magnetic field:', 'Hamiltonian model:']
        
        self.hamiltonian_edits = []
        for i, param in enumerate(hamiltonian_params[:-1]):
            label = QLabel(param)
            edit = QLineEdit()
            self.hamiltonian_edits.append(edit)
            hamiltonian_layout.addWidget(label, i, 0)
            hamiltonian_layout.addWidget(edit, i, 1)
        
        # Hamiltonian Model with Dropdown
        hamiltonian_label = QLabel(hamiltonian_params[-1])
        self.hamiltonian_dropdown = QComboBox()
        self.hamiltonian_dropdown.addItems(['Model 1', 'Model 2', 'Model 3'])
        hamiltonian_layout.addWidget(hamiltonian_label, len(hamiltonian_params) - 1, 0)
        hamiltonian_layout.addWidget(self.hamiltonian_dropdown, len(hamiltonian_params) - 1, 1)
        
        hamiltonian_group.setLayout(hamiltonian_layout)
        param_config_layout.addWidget(hamiltonian_group)

        # Qubit Topology Configuration (Subsection)
        qubit_group = QGroupBox('Qubit Topology Configuration')
        qubit_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 14px;
                font-weight: bold;
            }
        """)
        qubit_layout = QGridLayout()
        qubit_params = ['Lattice geometry:', 'Num. of spins:']
        
        # Lattice Geometry with Dropdown
        lattice_label = QLabel(qubit_params[0])
        self.lattice_dropdown = QComboBox()
        self.lattice_dropdown.addItems(['Geometry 1', 'Geometry 2', 'Geometry 3'])
        qubit_layout.addWidget(lattice_label, 0, 0)
        qubit_layout.addWidget(self.lattice_dropdown, 0, 1)
        
        # Number of Spins
        spins_label = QLabel(qubit_params[1])
        self.spins_edit = QLineEdit()
        qubit_layout.addWidget(spins_label, 1, 0)
        qubit_layout.addWidget(self.spins_edit, 1, 1)
        
        qubit_group.setLayout(qubit_layout)
        param_config_layout.addWidget(qubit_group)
        
        param_config_group.setLayout(param_config_layout)
        left_bottom_layout.addWidget(param_config_group)
        
        bottom_layout.addLayout(left_bottom_layout)
        
        # Right Section of Bottom Layout
        explore_group = QGroupBox('Explore')
        explore_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 16px;
                font-weight: bold;
            }
        """)
        explore_layout = QGridLayout()
        explore_params = ['Error difference:', 'Spin correlation:', 'Frustration level:']
        
        for i, param in enumerate(explore_params):
            label = QLabel(param)
            edit = QLineEdit()
            edit.setDisabled(True)  # Disable the input fields
            explore_layout.addWidget(label, i, 0)
            explore_layout.addWidget(edit, i, 1)
        
        explore_group.setLayout(explore_layout)
        bottom_layout.addWidget(explore_group)
        
        # Add bottom layout to main layout
        main_layout.addLayout(bottom_layout)

        # Submit Button
        submit_button = QPushButton('Submit', self)
        submit_button.setFixedSize(100, 40)  # Set size for the submit button
        submit_button.clicked.connect(self.onSubmit)
        param_config_layout.addWidget(submit_button, 0)  # Align center horizontally
        
        # Set the layout
        self.setLayout(main_layout)
        self.show()
        
    def onSubmit(self):
            # Gather input data
            interaction_strengths = self.hamiltonian_edits[0].text()
            magnetic_moments = self.hamiltonian_edits[1].text()
            external_magnetic_field = self.hamiltonian_edits[2].text()
            hamiltonian_model = self.hamiltonian_dropdown.currentText()
            lattice_geometry = self.lattice_dropdown.currentText()
            num_of_spins = int(self.spins_edit.text())
            
            # Create a quantum circuit based on user input
            qc = QuantumCircuit(num_of_spins)
            
            # Example: Apply some gates based on the Hamiltonian model
            if hamiltonian_model == 'Model 1':
                qc.h(0)
            elif hamiltonian_model == 'Model 2':
                qc.x(0)
            elif hamiltonian_model == 'Model 3':
                qc.y(0)
            
            # Draw the circuit
            circuit_image = circuit_drawer(qc, output='mpl')
            circuit_image.savefig('/tmp/circuit.png')  # Save the image temporarily

            # Load the image into the QLabel
            pixmap = QPixmap('/tmp/circuit.png')
            self.circuit_diagram_label.setPixmap(pixmap)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = QuantumCircuitGUI()
    sys.exit(app.exec())




