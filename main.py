import sys
import numpy as np
from matplotlib import pyplot as plt
from PyQt6.QtGui import QPixmap, QDoubleValidator
from qiskit import QuantumCircuit
from qiskit.visualization import circuit_drawer
from qiskit import *
from qiskit.circuit import Parameter
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit.primitives import BackendSampler
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QGroupBox, QGridLayout, QComboBox, QScrollArea, QCheckBox, QSlider
)

class QuantumCircuitGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.simulator = AerSimulator()
        self.initUI()
        self.Sx = []
        self.Sy = []
        self.Sz = []


    def initUI(self):
        """
        Initializes the user interface for the Quantum Circuit GUI.

        This function sets up the main window and its layout, including the circuit diagram group, scroll area, and circuit diagram labels. It also creates the parameter configuration group, spin system dropdown, and circuit configuration group. The function adds the necessary widgets and layouts to the main window and sets the initial visibility and default values of the spin system parameters.

        Parameters:
            None

        Returns:
            None
        """
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
        
        # Scroll area for circuit diagrams
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        
        # Widget to contain the circuit diagram labels
        self.circuit_widget = QWidget()
        circuit_diagram_layout = QVBoxLayout(self.circuit_widget)
        
        # Add QLabel widgets for the circuit diagrams
        self.circuit_diagram_label_1 = QLabel('', self)
        self.circuit_diagram_label_2 = QLabel('', self)
        self.circuit_diagram_label_3 = QLabel('', self)

        circuit_diagram_layout.addWidget(self.circuit_diagram_label_1)
        circuit_diagram_layout.addWidget(self.circuit_diagram_label_2)
        circuit_diagram_layout.addWidget(self.circuit_diagram_label_3)

        self.scroll_area.setWidget(self.circuit_widget)
        
        # Add the scroll area to the circuit diagram group layout
        circuit_diagram_group_layout = QVBoxLayout()
        circuit_diagram_group_layout.addWidget(self.scroll_area)
        self.circuit_diagram_group.setLayout(circuit_diagram_group_layout)

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

        spin_system_label = QLabel('Select the specific spin system:')
        self.spin_system_dropdown = QComboBox()
        self.spin_system_dropdown.addItems(['Spin-1', 'Spin-2', 'Spin-3'])
        self.spin_system_dropdown.currentIndexChanged.connect(self.update_spin_system_parameters)

        param_config_layout.addWidget(spin_system_label)
        param_config_layout.addWidget(self.spin_system_dropdown)

        # Hamiltonian Configuration (Subsection)
        hamiltonian_group = QGroupBox('Circuit Configuration')
        hamiltonian_group.setStyleSheet("""
            QGroupBox::title {
                font-size: 14px;
                font-weight: bold;
            }
        """)
        self.hamiltonian_layout = QGridLayout()
        circuit_params = ['Constant µ:', 
                        'Magnetic field strength:', 
                        'Time interval:', 
                        'Jx:', 
                        'Jy:', 
                        'Jz:', 
                        'theta_1:', 
                        'phi_1:', 
                        'lambda_1:', 
                        'theta_2:', 
                        'phi_2:', 
                        'lambda_2:', 
                        'theta_3:', 
                        'phi_3:', 
                        'lambda_3:', 
                        'Number of Trotter steps:']

        self.default_values = {
            'Spin-1': [1.0, 1.0, 0.01],  # Default values for Spin-1 system
            'Spin-2': [None, None, 0.1, 0.5, -0.45, 0.25, np.pi/6, np.pi/2, 0, 0, np.pi/4, np.pi/2],  # Default values for Spin-2 system
            'Spin-3': [None, None, 0.1, 0.5, -0.45, 0.25, np.pi/6, np.pi/2, 0, 0, np.pi/4, np.pi/2, 0, 0, 0, 100]  # Default values for Spin-3 system
        }

        self.hamiltonian_labels = []
        self.hamiltonian_edits = []
        self.double_validator = QDoubleValidator()
        for i, param in enumerate(circuit_params):
            label = QLabel(param)
            edit = QLineEdit()
            edit.setValidator(self.double_validator)  # Add validation
            self.hamiltonian_labels.append(label)
            self.hamiltonian_edits.append(edit)
            self.hamiltonian_layout.addWidget(label, i, 0)
            self.hamiltonian_layout.addWidget(edit, i, 1)

        hamiltonian_group.setLayout(self.hamiltonian_layout)
        param_config_layout.addWidget(hamiltonian_group)

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
        explore_layout = QVBoxLayout()

        # Slider
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(100)
        self.slider.setValue(50)
        self.slider.valueChanged.connect(self.on_slider_value_changed)
        explore_layout.addWidget(self.slider)
        self.time_step_label = QLabel('Time step: 0')
        explore_layout.addWidget(self.time_step_label)

        # Explore parameters
        explore_params = ['Sx:', 'Sy:', 'Sz:']

        self.explore_labels = []
        self.explore_edits = []

        for param in explore_params:
            label = QLabel(param)
            edit = QLineEdit()
            edit.setDisabled(True)  # Disable the input fields
            explore_layout.addWidget(label)
            explore_layout.addWidget(edit)
            self.explore_labels.append(label)
            self.explore_edits.append(edit)

        explore_group.setLayout(explore_layout)
        bottom_layout.addWidget(explore_group)

        # Add bottom layout to main layout
        main_layout.addLayout(bottom_layout)

        # Toggle switch and input field
        toggle_layout = QHBoxLayout()
        self.toggle_switch = QCheckBox('Enable API Key')
        self.api_key_input = QLineEdit()
        self.api_key_input.setPlaceholderText('Your API Key ...')
        self.api_key_input.setDisabled(True)

        self.toggle_switch.stateChanged.connect(self.on_toggle)

        toggle_layout.addWidget(self.toggle_switch)
        toggle_layout.addWidget(self.api_key_input)
        main_layout.addLayout(toggle_layout)

        # Submit Button
        submit_button = QPushButton('Submit', self)
        submit_button.setFixedSize(100, 40)  # Set size for the submit button
        submit_button.clicked.connect(self.onSubmit)
        param_config_layout.addWidget(submit_button, 0)  # Align center horizontally
        
        # Clear Button
        clear_button = QPushButton('Clear', self)
        clear_button.setFixedSize(100, 40)  # Set size for the clear button
        clear_button.clicked.connect(self.clear_circuits)
        param_config_layout.addWidget(clear_button, 0)  # Align center horizontally

        # Set the layout
        self.setLayout(main_layout)
        self.update_spin_system_parameters()  # Initialize the visibility of the parameters and set default values
        self.show()


    def build_circuit(self, params, selected_spin, is_api_key_enabled, api_key_value):
        """
        Builds the quantum circuit based on the selected spin system and parameters.

        Parameters:
            params (dict): A dictionary containing the parameters for the spin system.
            selected_spin (int): The index of the selected spin system (0 for Spin-1, 1 for Spin-2, 2 for Spin-3).
            is_api_key_enabled (bool): A flag indicating whether the API key is enabled.
            api_key_value (str): The value of the API key.

        Returns:
            None

        Description:
            This function builds the quantum circuit based on the selected spin system and parameters.
            It initializes the necessary variables and extracts the parameters from the `params` dictionary.
            Then, it calls the appropriate function to build the quantum circuit based on the selected spin system.
            Finally, it draws and saves the circuits, and updates the GUI with the circuit images.

        Note:
            - The function assumes that the necessary functions (`build_spin_one_system`, `build_spin_two_system`, `build_spin_three_system`, `plot_spin_one_system`, `plot_spin_two_system`, `plot_spin_three_system`) are defined elsewhere in the code.
            - The function assumes that the necessary GUI elements (`circuit_diagram_label_1`, `circuit_diagram_label_2`, `circuit_diagram_label_3`) are defined in the GUI class.
            - The function assumes that the necessary image files (`/tmp/circuit1.png`, `/tmp/circuit2.png`, `/tmp/circuit3.png`) are saved temporarily.

        """
        circ_1, circ_2, circ_3 = None, None, None
        spin_1, spin_2, spin_3 = False, False, False
        mu, B_0, del_t = 0, 0, 0
        Jx, Jy, Jz = 0, 0, 0
        theta_1, phi_1, lambda_1 = 0, 0, 0
        theta_2, phi_2, lambda_2 = 0, 0, 0
        theta_3, phi_3, lambda_3 = 0, 0, 0
        num_trotter_steps = 0

        if selected_spin == 0:  # If Spin-1 System is selected
            mu = float(params.get('Constant µ:', 0))
            B_0 = float(params.get('Magnetic field strength:', 1))
            del_t = float(params.get('Time interval:', 2))
            circ_1, (Sx, Sy, Sz) = self.build_spin_one_system(mu, B_0, del_t, is_api_key_enabled, api_key_value)
            self.Sx, self.Sy, self.Sz = Sx, Sy, Sz
            spin_1 = True

        elif selected_spin == 1:  # If Spin-2 System is selected
            del_t = float(params.get('Time interval:', 2))
            Jx = float(params.get('Jx:', 3))
            Jy = float(params.get('Jy:', 4))
            Jz = float(params.get('Jz:', 5))
            theta_1 = float(params.get('theta_1:', 6))
            phi_1 = float(params.get('phi_1:', 7))
            lambda_1 = float(params.get('lambda_1:', 8))
            theta_2 = float(params.get('theta_2:', 9))
            phi_2 = float(params.get('phi_2:', 10))
            lambda_2 = float(params.get('lambda_2:', 11))

            circ_1, circ_2, circ_3, (Sx1, Sy1), (Sz1, Sx2), (Sy2, Sz2) = self.build_spin_two_system(del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, is_api_key_enabled, api_key_value)
            self.Sx, self.Sy, self.Sz = Sx1, Sy1, Sz1
            spin_2 = True

        elif selected_spin == 2:  # If Spin-3 System is selected
            del_t = float(params.get('Time interval:', 2))
            Jx = float(params.get('Jx:', 3))
            Jy = float(params.get('Jy:', 4))
            Jz = float(params.get('Jz:', 5))
            theta_1 = float(params.get('theta_1:', 6))
            phi_1 = float(params.get('phi_1:', 7))
            lambda_1 = float(params.get('lambda_1:', 8))
            theta_2 = float(params.get('theta_2:', 9))
            phi_2 = float(params.get('phi_2:', 10))
            lambda_2 = float(params.get('lambda_2:', 11))
            theta_3 = float(params.get('theta_3:', 12))
            phi_3 = float(params.get('phi_3:', 13))
            lambda_3 = float(params.get('lambda_3:', 14))
            num_trotter_steps = int(params.get('Number of Trotter steps:', 15))

            circ_1, circ_2, circ_3, Sx, Sy, Sz = self.build_spin_three_system(del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, num_trotter_steps, is_api_key_enabled, api_key_value)
            self.Sx, self.Sy, self.Sz = Sx, Sy, Sz
            spin_3 = True

        # Draw and save the circuits
        if spin_1:
            Nt = int(2*np.pi // (mu*B_0*del_t))
            tau_range = np.linspace(0, 2*np.pi, Nt)

            self.plot_spin_one_system(tau_range, Sx, Sy, Sz, Nt)

            circuit_image_1 = circuit_drawer(circ_1, output='mpl')
            circuit_image_1.savefig('/tmp/circuit1.png')  # Save the image temporarily
            pixmap_1 = QPixmap('/tmp/circuit1.png')
            self.circuit_diagram_label_1.setPixmap(pixmap_1)

            spin_1 = False

        if spin_2:
            Nt = int(10//del_t)
            tau_range = np.linspace(0,10,Nt)

            self.plot_spin_two_system(tau_range, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2)

            circuit_image_1 = circuit_drawer(circ_1, output='mpl')
            circuit_image_1.savefig('/tmp/circuit1.png')  # Save the image temporarily
            pixmap_1 = QPixmap('/tmp/circuit1.png')
            self.circuit_diagram_label_1.setPixmap(pixmap_1)

            circuit_image_2 = circuit_drawer(circ_2, output='mpl')
            circuit_image_2.savefig('/tmp/circuit2.png')  # Save the image temporarily
            pixmap_2 = QPixmap('/tmp/circuit2.png')
            self.circuit_diagram_label_2.setPixmap(pixmap_2)

            circuit_image_3 = circuit_drawer(circ_3, output='mpl')
            circuit_image_3.savefig('/tmp/circuit3.png')  # Save the image temporarily
            pixmap_3 = QPixmap('/tmp/circuit3.png')
            self.circuit_diagram_label_3.setPixmap(pixmap_3)

            spin_2 = False

        if spin_3:
            Nt = int(10//del_t)
            tau_range = np.linspace(0,10,Nt)

            self.plot_spin_three_system(tau_range, Sx, Sy, Sz)

            circuit_image_1 = circuit_drawer(circ_1, output='mpl')
            circuit_image_1.savefig('/tmp/circuit1.png')  # Save the image temporarily
            pixmap_1 = QPixmap('/tmp/circuit1.png')
            self.circuit_diagram_label_1.setPixmap(pixmap_1)

            circuit_image_2 = circuit_drawer(circ_2, output='mpl')
            circuit_image_2.savefig('/tmp/circuit2.png')  # Save the image temporarily
            pixmap_2 = QPixmap('/tmp/circuit2.png')
            self.circuit_diagram_label_2.setPixmap(pixmap_2)

            circuit_image_3 = circuit_drawer(circ_3, output='mpl')
            circuit_image_3.savefig('/tmp/circuit3.png')  # Save the image temporarily
            pixmap_3 = QPixmap('/tmp/circuit3.png')
            self.circuit_diagram_label_3.setPixmap(pixmap_3)

            spin_3 = False   
    def build_spin_one_system(self, mu, B_0, del_t, is_api_key_enabled, api_key_value):
        """
        Builds a quantum circuit for simulating a spin one system.

        Args:
            mu (float): The magnetic field magnitude.
            B_0 (float): The static magnetic field.
            del_t (float): The time step.
            is_api_key_enabled (bool): Whether the API key is enabled.
            api_key_value (str): The value of the API key.

        Returns:
            tuple: A tuple containing the built quantum circuit and the result of running the circuit.

        Raises:
            None
        """
        tau = Parameter('τ')
        qr = QuantumRegister(3,'q')
        cr = ClassicalRegister(3,'c')

        Nt = int(2*np.pi // (mu*B_0*del_t))

        tau_range = np.linspace(0, 2*np.pi, Nt)

        timecirc = QuantumCircuit(qr,cr)

        #no initial unitary transformation since |ψ(0)> = |+>

        timecirc.u(tau,np.pi/2,-np.pi/2,qr) #apply exp(-iHt/ħ)
        timecirc.barrier(qr)
        timecirc.ry(-np.pi/2,0) #rotation to measure <Sx>

        timecirc.rz(-np.pi/2,1)
        timecirc.ry(-np.pi/2,1) #rotation to measure <Sy>
        timecirc.barrier(qr)
        #no rotation needed to measure <Sz>

        timecirc.measure(qr,cr)

        return timecirc, self.run_spin_one_system(timecirc, tau_range, tau, self.simulator, is_api_key_enabled, api_key_value)
    def build_spin_two_system(self, del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, is_api_key_enabled, api_key_value):
        """
        Builds and returns three quantum circuits for simulating the spin two-system problem.
        
        Args:
            del_t (float): The time interval.
            Jx (float): The x-axis coupling strength.
            Jy (float): The y-axis coupling strength.
            Jz (float): The z-axis coupling strength.
            theta_1 (float): The angle for the first spin system.
            phi_1 (float): The phase for the first spin system.
            lambda_1 (float): The phase for the first spin system.
            theta_2 (float): The angle for the second spin system.
            phi_2 (float): The phase for the second spin system.
            lambda_2 (float): The phase for the second spin system.
            is_api_key_enabled (bool): Flag indicating whether the API key is enabled.
            api_key_value (str): The API key value.
            
        Returns:
            tuple: A tuple containing three QuantumCircuit objects representing the Sx, Sy, and Sz circuits.
            tuple: A tuple containing the results of running the Sx, Sy, and Sz circuits using the specified time range and parameter.
        """
        Nt = int(10//del_t)
        tau = Parameter('τ')
        tau_range = np.linspace(0,10,Nt)
        
        def build_Sx():
            qr = QuantumRegister(2)
            cr = ClassicalRegister(2)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)

            self.N((Jx*tau/4.0),(Jy*tau/4.0),(Jz*tau/4.0),timecirc,0,1)

            timecirc.ry(-np.pi/2,qr)
            timecirc.measure(qr,[1,0])
            
            return timecirc
        def build_Sy():
            qr = QuantumRegister(2)
            cr = ClassicalRegister(2)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)

            self.N(Jx*tau/4.0,Jy*tau/4.0,Jz*tau/4.0,timecirc,0,1)

            timecirc.rz(-np.pi/2,qr)
            timecirc.ry(-np.pi/2,qr)
            timecirc.measure(qr,[1,0])

            return timecirc
        def build_Sz():
            qr = QuantumRegister(2)
            cr = ClassicalRegister(2)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)

            self.N(Jx*tau/4.0,Jy*tau/4.0,Jz*tau/4.0,timecirc,0,1)

            timecirc.measure(qr,[1,0])

            return timecirc
        
        Sx_circ = build_Sx()
        Sy_circ = build_Sy()
        Sz_circ = build_Sz()
        
        return Sx_circ, Sy_circ, Sz_circ, self.run_spin_two_system(Sx_circ, tau_range, tau, self.simulator, is_api_key_enabled, api_key_value), self.run_spin_two_system(Sy_circ, tau_range, tau, self.simulator, is_api_key_enabled, api_key_value), self.run_spin_two_system(Sz_circ, tau_range, tau, self.simulator, is_api_key_enabled, api_key_value)
    def build_spin_three_system(self, del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, ntrot, is_api_key_enabled, api_key_value):
        """
        Builds the spin three system quantum circuit for measuring the spin operators Sx, Sy, and Sz.
        
        Args:
            del_t (float): The time step for the simulation.
            Jx (float): The x-axis exchange interaction strength.
            Jy (float): The y-axis exchange interaction strength.
            Jz (float): The z-axis exchange interaction strength.
            theta_1 (float): The theta parameter for the first spin state.
            phi_1 (float): The phi parameter for the first spin state.
            lambda_1 (float): The lambda parameter for the first spin state.
            theta_2 (float): The theta parameter for the second spin state.
            phi_2 (float): The phi parameter for the second spin state.
            lambda_2 (float): The lambda parameter for the second spin state.
            theta_3 (float): The theta parameter for the third spin state.
            phi_3 (float): The phi parameter for the third spin state.
            lambda_3 (float): The lambda parameter for the third spin state.
            ntrot (int): The number of Trotter steps.
            is_api_key_enabled (bool): Indicates whether the API key is enabled.
            api_key_value (str): The value of the API key.
        
        Returns:
            tuple: A tuple containing the quantum circuits for measuring Sx, Sy, and Sz, and the simulation results for each operator.
                - Sx_circ (QuantumCircuit): The quantum circuit for measuring Sx.
                - Sy_circ (QuantumCircuit): The quantum circuit for measuring Sy.
                - Sz_circ (QuantumCircuit): The quantum circuit for measuring Sz.
                - Sx_result (list): The simulation results for Sx.
                - Sy_result (list): The simulation results for Sy.
                - Sz_result (list): The simulation results for Sz.
        """
        Nt = int(10//del_t)
        tau = Parameter('τ')
        tau_range = np.linspace(0,10,Nt)

        def build_Sx():
            qr = QuantumRegister(3)
            cr = ClassicalRegister(3)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)
            timecirc.u(theta_3,phi_3,lambda_3,2)

            for i in range(0,ntrot):
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,0,1)
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,1,2)

            timecirc.ry(-np.pi/2,qr)

            timecirc.measure(qr,cr)
            
            return timecirc
        def build_Sy():
            qr = QuantumRegister(3)
            cr = ClassicalRegister(3)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)
            timecirc.u(theta_3,phi_3,lambda_3,2)

            for i in range(0,ntrot):
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,0,1)
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,1,2)

            timecirc.rz(-np.pi/2,qr)
            timecirc.ry(-np.pi/2,qr)
            timecirc.measure(qr,cr)

            return timecirc
        def build_Sz():
            qr = QuantumRegister(3)
            cr = ClassicalRegister(3)

            timecirc = QuantumCircuit(qr,cr)

            #initial states
            timecirc.u(theta_1,phi_1,lambda_1,0)
            timecirc.u(theta_2,phi_2,lambda_2,1)
            timecirc.u(theta_3,phi_3,lambda_3,2)

            for i in range(0,ntrot):
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,0,1)
                self.N(Jx*(tau/ntrot)/4.0,Jy*(tau/ntrot)/4.0,Jz*(tau/ntrot)/4.0,timecirc,1,2)

            timecirc.measure(qr,cr)
            
            return timecirc

        Sx_circ = build_Sx()
        Sy_circ = build_Sy()
        Sz_circ = build_Sz()

        return Sx_circ, Sy_circ, Sz_circ, self.run_spin_three_system(Sx_circ, tau_range, tau, self.simulator, Nt, is_api_key_enabled, api_key_value), self.run_spin_three_system(Sy_circ, tau_range, tau, self.simulator, Nt, is_api_key_enabled, api_key_value), self.run_spin_three_system(Sz_circ, tau_range, tau, self.simulator, Nt, is_api_key_enabled, api_key_value)

    def run_spin_one_system(self, timecirc, tau_range, tau, simulator, is_api_key_enabled, api_key_value):
        """
        Runs the spin one system simulation using the provided quantum circuit, tau range, tau parameter, simulator, API key status, and API key value.

        Parameters:
            timecirc (QuantumCircuit): The quantum circuit to be simulated.
            tau_range (list): The range of tau values to iterate over.
            tau (Parameter): The tau parameter to be assigned to the quantum circuit.
            simulator (Simulator): The simulator to be used for simulation.
            is_api_key_enabled (bool): Indicates whether the API key is enabled.
            api_key_value (str): The value of the API key.

        Returns:
            list: The postprocessed spin one system simulation results.
        """
        simcounts = []
        if is_api_key_enabled is False and api_key_value == "":
            Nshots = 8192
            transpiled_circ = transpile(timecirc, simulator)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()
                simcounts.append(result.get_counts(transpiled_circ_with_param))
        else:
            service, backend = self.setup_quantum_hardware(api_key_value)
            sampler = BackendSampler(backend=backend)
            transpiled_circ = transpile(timecirc, backend)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                job_with_result = sampler.run([transpiled_circ_with_param]).result()
                simcounts.append(job_with_result.quasi_dists[0])

        return self.postprocess_spin_one_system(simcounts, Nshots)
    def run_spin_two_system(self, timecirc, tau_range, tau, simulator, is_api_key_enabled, api_key_value):
        """
        Runs a quantum circuit simulation for a two-system spin system.

        Args:
            timecirc (QuantumCircuit): The quantum circuit to simulate.
            tau_range (list): The range of values for the tau parameter.
            tau (Parameter): The tau parameter.
            simulator (Simulator): The simulator to use for simulation.
            is_api_key_enabled (bool): Flag indicating whether the API key is enabled.
            api_key_value (str): The value of the API key.

        Returns:
            list: The postprocessed simulation counts.

        Raises:
            None
        """
        simcounts = []
        if is_api_key_enabled is False and api_key_value == "":
            Nshots = 8192
            transpiled_circ = transpile(timecirc, simulator)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()
                simcounts.append(result.get_counts(transpiled_circ_with_param))
        else:
            service, backend = self.setup_quantum_hardware(api_key_value)
            sampler = BackendSampler(backend=backend)
            transpiled_circ = transpile(timecirc, backend)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                job_with_result = sampler.run([transpiled_circ_with_param]).result()
                simcounts.append(job_with_result.quasi_dists[0])

        return self.postprocess_spin_two_system(simcounts, Nshots)
    def run_spin_three_system(self, timecirc, tau_range, tau, simulator, Nt, is_api_key_enabled, api_key_value):
        """
        Runs a quantum circuit simulation for a three-system spin system.

        Args:
            timecirc (QuantumCircuit): The quantum circuit to be simulated.
            tau_range (Iterable): The range of tau values to iterate over.
            tau (Parameter): The parameter representing tau in the quantum circuit.
            simulator (Simulator): The simulator to use for the simulation.
            Nt (int): The number of time steps in the simulation.
            is_api_key_enabled (bool): Flag indicating whether the API key is enabled.
            api_key_value (str): The value of the API key.

        Returns:
            list: The postprocessed simulation counts for each tau value.

        Description:
            This function runs a quantum circuit simulation for a three-system spin system. It takes in a quantum circuit,
            a range of tau values, a parameter representing tau, a simulator, the number of time steps in the simulation,
            a flag indicating whether the API key is enabled, and the value of the API key. It then performs the simulation
            for each tau value in the range. If the API key is not enabled and the API key value is empty, it uses the provided
            simulator to run the simulation with 8192 shots. Otherwise, it sets up the quantum hardware using the API key value
            and uses the BackendSampler to run the simulation. The simulation counts for each tau value are stored in a list
            and returned after postprocessing.
        """
        simcounts = []
        if is_api_key_enabled is False and api_key_value == "":
            Nshots = 8192
            transpiled_circ = transpile(timecirc, simulator)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()
                simcounts.append(result.get_counts(transpiled_circ_with_param))
        else:
            service, backend = self.setup_quantum_hardware(api_key_value)
            sampler = BackendSampler(backend=backend)
            transpiled_circ = transpile(timecirc, backend)
            for t in tau_range:
                transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})
                job_with_result = sampler.run([transpiled_circ_with_param]).result()
                simcounts.append(job_with_result.quasi_dists[0])

        return self.postprocess_spin_three_system(simcounts, Nshots, Nt)

    def postprocess_spin_one_system(self, simcounts, Nshots):
        """
        Calculates the spin operators Sx, Sy, and Sz for a one-system spin chain.

        Parameters:
            simcounts (list): A list of dictionaries containing the simulation counts for each qubit.
            Nshots (int): The number of shots taken in the simulation.

        Returns:
            tuple: A tuple containing the values of Sx, Sy, and Sz.
        """
        c000 = np.array(list(map(lambda c: c.get('000', 0), simcounts)))
        c001 = np.array(list(map(lambda c: c.get('001', 0), simcounts)))
        c010 = np.array(list(map(lambda c: c.get('010', 0), simcounts)))
        c011 = np.array(list(map(lambda c: c.get('011', 0), simcounts)))
        c100 = np.array(list(map(lambda c: c.get('100', 0), simcounts)))
        c101 = np.array(list(map(lambda c: c.get('101', 0), simcounts)))
        c110 = np.array(list(map(lambda c: c.get('110', 0), simcounts)))
        c111 = np.array(list(map(lambda c: c.get('111', 0), simcounts)))

        Sz = 0.5*(c000+c001+c010+c011-c100-c101-c110-c111)/Nshots
        Sy = 0.5*(c000+c001+c100+c101-c010-c011-c110-c111)/Nshots
        Sx = 0.5*(c000+c010+c100+c110-c001-c011-c101-c111)/Nshots

        return Sx, Sy, Sz
    def postprocess_spin_two_system(self, simcounts, Nshots):
        """
        Calculates the spin operators S1 and S2 for a two-system spin chain.

        Parameters:
            simcounts (list): A list of dictionaries containing the simulation counts for each qubit.
            Nshots (int): The number of shots taken in the simulation.

        Returns:
            tuple: A tuple containing the values of S1 and S2.
        """
        c00 = np.array(list(map(lambda c: c.get('00', 0), simcounts)))
        c01 = np.array(list(map(lambda c: c.get('01', 0), simcounts)))
        c10 = np.array(list(map(lambda c: c.get('10', 0), simcounts)))
        c11 = np.array(list(map(lambda c: c.get('11', 0), simcounts)))

        S1 = 0.5*(c00+c01-c10-c11)/Nshots
        S2 = 0.5*(c00+c10-c01-c11)/Nshots

        return S1, S2 
    def postprocess_spin_three_system(self, simcounts, Nshots, Nt):
        """
        Postprocesses the simulation counts for a three-system spin chain.

        Args:
            simcounts (list): A list of dictionaries containing the simulation counts for each qubit.
            Nshots (int): The number of shots taken in the simulation.
            Nt (int): The number of time steps in the simulation.

        Returns:
            numpy.ndarray: An array of shape (3, Nt) containing the time-dependent expectation values of the spin operators S1, S2, and S3.
        """
        St = np.zeros((3,Nt))
        for i in range(0,Nt):
            counts = simcounts[i]
            keylist = list(counts.keys())
            for j in range(0,len(keylist)):
                state = keylist[j]

                if (state[0]=='0'):
                    St[0,i] = St[0,i] + counts[state]
                else:
                    St[0,i] = St[0,i] - counts[state]


                if (state[1]=='0'):
                    St[1,i] = St[1,i] + counts[state]
                else:
                    St[1,i] = St[1,i] - counts[state]


                if (state[2]=='0'):
                    St[2,i] = St[2,i] + counts[state]
                else:
                    St[2,i] = St[2,i] - counts[state]

        St = 0.5*St/Nshots

        return St

    def plot_spin_one_system(self, tau_range, Sx, Sy, Sz, Nt):
        """
        Plots the time evolution graph of a single spin.

        Args:
            tau_range (array-like): The range of values for the tau parameter.
            Sx (array-like): The values of the Sx expectation value.
            Sy (array-like): The values of the Sy expectation value.
            Sz (array-like): The values of the Sz expectation value.
            Nt (int): The number of time steps.

        Returns:
            None

        This function creates a new figure and plots the time evolution graph of a single spin. It takes in the range of values for the tau parameter, the values of the Sx, Sy, and Sz expectation values, and the number of time steps. The function then plots the Sx, Sy, and Sz expectation values against the tau parameter using different colors and styles. It also plots a constant line and a sine and cosine wave for comparison. The x-axis is labeled as "ωt" and the y-axis is labeled as "<S>". The function also adds a legend with labels for the different plots. Finally, it adds a title to the plot and displays it using plt.show().
        """
        plt.close('all')  # Close all open figure windows
        plt.figure()  # Create a new figure
        plt.plot(tau_range, Sx, 'bo', label='<Sx>')
        plt.plot(tau_range, Sy, 'ro', label='<Sy>')
        plt.plot(tau_range, Sz, 'ko', label='<Sz>')
        plt.plot(tau_range, 0 * np.zeros(Nt), 'b-')
        plt.plot(tau_range, 0.5 * np.sin(tau_range), 'r-')
        plt.plot(tau_range, 0.5 * np.cos(tau_range), 'k-')
        plt.xlabel('$\omega t$')
        plt.ylabel('<S>')
        plt.legend()
        plt.title("Time evolution graph of a single spin")
        plt.show()
    def plot_spin_two_system(self, tau_range, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2):
        """
        Plots the time evolution graph of two spin systems.

        Args:
            tau_range (array-like): The range of tau values.
            Sx1 (array-like): The values of Sx1 at each tau value.
            Sy1 (array-like): The values of Sy1 at each tau value.
            Sz1 (array-like): The values of Sz1 at each tau value.
            Sx2 (array-like): The values of Sx2 at each tau value.
            Sy2 (array-like): The values of Sy2 at each tau value.
            Sz2 (array-like): The values of Sz2 at each tau value.

        Returns:
            None

        This function plots the time evolution graph of two spin systems. It creates two figures, one for each spin system. The x-axis represents the tau values, and the y-axis represents the values of Sx, Sy, and Sz at each tau value. The legend displays the labels 'j=x', 'j=y', and 'j=z' for each spin system. The title of each figure provides a brief description of the spin system being plotted.

        Note:
            - The function uses the matplotlib library to create the plots.
            - The function assumes that the input arrays have the same length.
            - The function does not return any values.
        """
        # Time evolution graph of the first spin
        plt.close('all')  # Close all open figure windows
        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sx1,'r--')
        plt.plot(tau_range,Sy1,'b-.')
        plt.plot(tau_range,Sz1,'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.5,0.5)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.title("Time evolution graph of the first spin")
        plt.show()

        # Time evolution graph of the second spin
        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sx2,'r--')
        plt.plot(tau_range,Sy2,'b-.')
        plt.plot(tau_range,Sz2,'k-')
        plt.xlabel('Jt')
        plt.ylabel('$<S_{2}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.title("Time evolution graph of the second spin")
        plt.show()
    def plot_spin_three_system(self, tau_range, Sxt, Syt, Szt):
        """
        Plots the time evolution graph of three spin systems.

        Args:
            tau_range (array-like): The range of tau values.
            Sxt (array-like): The values of Sxt at each tau value for each spin system.
            Syt (array-like): The values of Syt at each tau value for each spin system.
            Szt (array-like): The values of Szt at each tau value for each spin system.

        Returns:
            None

        This function plots the time evolution graph of three spin systems. It creates three figures, one for each spin system. The x-axis represents the tau values, and the y-axis represents the values of Sxt, Syt, and Szt at each tau value. The legend displays the labels 'j=x', 'j=y', and 'j=z' for each spin system. The title of each figure provides a brief description of the spin system being plotted.

        Note:
            - The function uses the matplotlib library to create the plots.
            - The function assumes that the input arrays have the same length.
            - The function does not return any values.
        """
        # Time evolution graph of the first spin
        plt.close('all')  # Close all open figure windows
        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sxt[0,:],'r--')
        plt.plot(tau_range,Syt[0,:],'b-.')
        plt.plot(tau_range,Szt[0,:],'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.6,0.6)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.title("Time evolution graph of the first spin")
        plt.show()

        # Time evolution graph of the second spin
        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sxt[1,:],'r--')
        plt.plot(tau_range,Syt[1,:],'b-.')
        plt.plot(tau_range,Szt[1,:],'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.6,0.6)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.title("Time evolution graph of the second spin")
        plt.show()

        # Time evolution graph of the third spin
        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sxt[2,:],'r--')
        plt.plot(tau_range,Syt[2,:],'b-.')
        plt.plot(tau_range,Szt[2,:],'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.6,0.6)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.title("Time evolution graph of the third spin")
        plt.show()


    def N(self, α, β, γ, circ, q1, q2):
        """
        Applies a series of quantum gates to a given quantum circuit to simulate the action of the N operator.

        Parameters:
            α (float): The angle parameter for the first rotation gate.
            β (float): The angle parameter for the second rotation gate.
            γ (float): The angle parameter for the third rotation gate.
            circ (QuantumCircuit): The quantum circuit to apply the gates to.
            q1 (int): The index of the first qubit.
            q2 (int): The index of the second qubit.

        Returns:
            None
        """
        circ.rz(-0.5*np.pi,q2)
        circ.cx(q2,q1)
        circ.rz(0.5*np.pi-2*γ,q1)
        circ.ry(2.0*α-0.5*np.pi,q2)
        circ.cx(q1,q2)
        circ.ry(0.5*np.pi-2.0*β,q2)
        circ.cx(q2,q1)
        circ.rz(0.5*np.pi,q1)
    def clear_circuits(self):
        """
        Clears the circuit diagrams displayed in the GUI by clearing the corresponding QLabel widgets.

        This function does not take any parameters.

        This function does not return anything.
        """
        self.circuit_diagram_label_1.clear()
        self.circuit_diagram_label_2.clear()
        self.circuit_diagram_label_3.clear()
    def update_spin_system_parameters(self):
        """
        Updates the parameters of the spin system based on the selected spin type.

        This function retrieves the selected spin type from the `spin_system_dropdown` widget and 
        determines the corresponding parameters to be shown. It then iterates over the labels and 
        edits of the `hamiltonian_labels` and `hamiltonian_edits` lists respectively, and shows or 
        hides them based on whether they are in the list of parameters to be shown. If the parameter 
        is to be shown, it also sets the initial value of the edit widget using the corresponding 
        default value from the `self.default_values` dictionary.

        After updating the visibility and values of the labels and edits, the function disables the 
        `slider` widget if the selected spin type is Spin-3, otherwise it enables it.

        Parameters:
            None

        Returns:
            None
        """
        self.selected_spin = self.spin_system_dropdown.currentIndex()
        params_to_show = {
            0: [0, 1, 2],  # Spin-1: Constant µ, Magnetic field strength, Time interval
            1: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],  # Spin-2: Time interval, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2
            2: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]  # Spin-3: Time interval, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, Number of Trotter steps
        }
        
        spin_types = ['Spin-1', 'Spin-2', 'Spin-3']
        spin_type = spin_types[self.selected_spin]
        
        for i, (label, edit) in enumerate(zip(self.hamiltonian_labels, self.hamiltonian_edits)):
            if i in params_to_show[self.selected_spin]:
                label.show()
                edit.show()
                # Set initial values
                if i < len(self.default_values[spin_type]):
                    edit.setText(str(self.default_values[spin_type][i]))
            else:
                label.hide()
                edit.hide()

            # Disable the slider if Spin-3 is selected, otherwise enable it
        if self.selected_spin == 2:
            self.slider.setDisabled(True)
        else:
            self.slider.setDisabled(False)
    def onSubmit(self):
        """
        Gathers input data and calls the `build_circuit` method with the gathered data.

        This function retrieves the input data from the UI and prepares it for the `build_circuit` method. It first creates an empty dictionary `circuit_params` to store the gathered input data. Then, it iterates over the `self.hamiltonian_edits` list and checks if each edit is visible. If it is, it adds the corresponding label text as the key and the edit text as the value to the `circuit_params` dictionary.

        After gathering the input data, the function checks the state of the `self.toggle_switch` widget to determine if the API key is enabled. It retrieves the value of the `self.api_key_input` widget.

        Finally, the function calls the `build_circuit` method with the gathered `circuit_params`, `self.selected_spin`, `is_api_key_enabled`, and `api_key_value` as arguments.

        Parameters:
            self (object): The instance of the class.
        
        Returns:
            None
        """
        # Gather input data
        circuit_params = {}
        for i, edit in enumerate(self.hamiltonian_edits):
            if edit.isVisible():
                circuit_params[self.hamiltonian_labels[i].text()] = edit.text()
        
        is_api_key_enabled = self.toggle_switch.isChecked()
        api_key_value = self.api_key_input.text()

        self.build_circuit(circuit_params, self.selected_spin, is_api_key_enabled, api_key_value)
    def setup_quantum_hardware(self, token):
        """
        Sets up the quantum hardware for use with Qiskit.

        Args:
            token (str): The authentication token for accessing the IBM Quantum API.

        Returns:
            tuple: A tuple containing the QiskitRuntimeService object and the least busy backend that meets the specified criteria.

        Raises:
            None

        Example Usage:
            service, backend = setup_quantum_hardware(token)
        """
        service = QiskitRuntimeService(channel="ibm_quantum", token=token)
        backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
        print("Backend online!")
        return service, backend
    def on_toggle(self):
        """
        Toggles the state of the API key input field based on the state of the toggle switch.

        This function is called when the toggle switch is clicked. It checks the state of the toggle switch
        and enables or disables the API key input field accordingly. If the toggle switch is checked, the
        API key input field is enabled, allowing the user to enter an API key. If the toggle switch is not
        checked, the API key input field is disabled, preventing the user from entering an API key.

        Parameters:
            self (object): The instance of the class.

        Returns:
            None
        """
        if self.toggle_switch.isChecked():
            self.api_key_input.setDisabled(False)
        else:
            self.api_key_input.setDisabled(True)
    def on_slider_value_changed(self):
        """
        Handles the event when the slider value is changed.

        This function is called when the slider value is changed by the user. It retrieves the current slider value, calculates the corresponding index based on the length of the Sx list, and updates the text of the explore edits and time step label accordingly.

        Parameters:
            self (object): The instance of the class.
        
        Returns:
            None
        """
        slider_value = self.slider.value()
        Nt = len(self.Sx)
        if Nt > 0:
            index = int((slider_value / 100) * (Nt - 1))
            print(self.Sx[index])
            self.explore_edits[0].setText(f'{self.Sx[index]:.4f}')
            self.explore_edits[1].setText(f'{self.Sy[index]:.4f}')
            self.explore_edits[2].setText(f'{self.Sz[index]:.4f}')
            self.time_step_label.setText(f'Time step: {index}')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = QuantumCircuitGUI()
    sys.exit(app.exec())
