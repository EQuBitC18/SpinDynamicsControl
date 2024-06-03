import os
import sys
import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QLineEdit, QPushButton, QGroupBox, QGridLayout, QComboBox, QScrollArea
)
from PyQt6.QtGui import QPixmap, QDoubleValidator
from qiskit import QuantumCircuit
from qiskit.visualization import circuit_drawer
from matplotlib import pyplot as plt
from qiskit import *
from qiskit.circuit import Parameter
from qiskit_aer import AerSimulator

class QuantumCircuitGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.simulator = AerSimulator()
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
        
        # Clear Button
        clear_button = QPushButton('Clear', self)
        clear_button.setFixedSize(100, 40)  # Set size for the clear button
        clear_button.clicked.connect(self.clear_circuits)
        param_config_layout.addWidget(clear_button, 0)  # Align center horizontally

        # Set the layout
        self.setLayout(main_layout)
        self.update_spin_system_parameters()  # Initialize the visibility of the parameters
        self.show()


    def build_circuit(self, params, selected_spin):
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
            circ_1, (Sx, Sy, Sz) = self.build_spin_one_system(mu, B_0, del_t)
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

            circ_1, circ_2, circ_3, (Sx1, Sy1), (Sz1, Sx2), (Sy2, Sz2) = self.build_spin_two_system(del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2)
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

            circ_1, circ_2, circ_3, Sx, Sy, Sz = self.build_spin_three_system(del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, num_trotter_steps)
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
            circ_1 = None, None, None

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

            circuit_image_3 = circuit_drawer(circ_3, output='mpl')
            circuit_image_3.savefig('/tmp/circuit3.png')  # Save the image temporarily
            pixmap_3 = QPixmap('/tmp/circuit3.png')
            self.circuit_diagram_label_3.setPixmap(pixmap_3)

            spin_3 = False   
    def build_spin_one_system(self, mu, B_0, del_t):
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

        return timecirc, self.run_spin_one_system(timecirc, tau_range, tau, self.simulator)
    def build_spin_two_system(self, del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2):
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
        
        return Sx_circ, Sy_circ, Sz_circ, self.run_spin_two_system(Sx_circ, tau_range, tau, self.simulator), self.run_spin_two_system(Sy_circ, tau_range, tau, self.simulator), self.run_spin_two_system(Sz_circ, tau_range, tau, self.simulator)
    def build_spin_three_system(self, del_t, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, ntrot):
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

        return Sx_circ, Sy_circ, Sz_circ, self.run_spin_three_system(Sx_circ, tau_range, tau, self.simulator, Nt), self.run_spin_three_system(Sy_circ, tau_range, tau, self.simulator, Nt), self.run_spin_three_system(Sz_circ, tau_range, tau, self.simulator, Nt)

    def run_spin_one_system(self, timecirc, tau_range, tau, simulator):
        Nshots = 8192
        simcounts = []
        transpiled_circ = transpile(timecirc, simulator)
        for t in tau_range:
            transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})

            result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()

            simcounts.append(result.get_counts(transpiled_circ_with_param))

        return self.postprocess_spin_one_system(simcounts, Nshots)
    def run_spin_two_system(self, timecirc, tau_range, tau, simulator):
        transpiled_circ = transpile(timecirc, simulator)

        Nshots = 8192
        simcounts = []
        for t in tau_range:
            transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})

            result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()

            simcounts.append(result.get_counts(transpiled_circ_with_param))

        return self.postprocess_spin_two_system(simcounts, Nshots)
    def run_spin_three_system(self, timecirc, tau_range, tau, simulator, Nt):
        transpiled_circ = transpile(timecirc, simulator)

        Nshots = 8192
        simcounts = []
        for t in tau_range:
            transpiled_circ_with_param = transpiled_circ.assign_parameters({tau: t})

            result = simulator.run(transpiled_circ_with_param, shots=Nshots).result()

            simcounts.append(result.get_counts(transpiled_circ_with_param))

        return self.postprocess_spin_three_system(simcounts, Nshots, Nt)

    def postprocess_spin_one_system(self, simcounts, Nshots):
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
        c00 = np.array(list(map(lambda c: c.get('00', 0), simcounts)))
        c01 = np.array(list(map(lambda c: c.get('01', 0), simcounts)))
        c10 = np.array(list(map(lambda c: c.get('10', 0), simcounts)))
        c11 = np.array(list(map(lambda c: c.get('11', 0), simcounts)))

        S1 = 0.5*(c00+c01-c10-c11)/Nshots
        S2 = 0.5*(c00+c10-c01-c11)/Nshots

        return S1, S2 
    def postprocess_spin_three_system(self, simcounts, Nshots, Nt):
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
        plt.show()
    def plot_spin_two_system(self, tau_range, Sx1, Sy1, Sz1, Sx2, Sy2, Sz2):
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
        plt.show()

        plt.figure()  # Create a new figure
        plt.plot(tau_range,Sx2,'r--')
        plt.plot(tau_range,Sy2,'b-.')
        plt.plot(tau_range,Sz2,'k-')
        plt.xlabel('Jt')
        plt.ylabel('$<S_{2}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.show()
    def plot_spin_three_system(self, tau_range, Sxt, Syt, Szt):
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
        plt.show()

        # Time evolution graph of the second spin
        plt.plot(tau_range,Sxt[1,:],'r--')
        plt.plot(tau_range,Syt[1,:],'b-.')
        plt.plot(tau_range,Szt[1,:],'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.6,0.6)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.show()

        # Time evolution graph of the third spin
        plt.plot(tau_range,Sxt[2,:],'r--')
        plt.plot(tau_range,Syt[2,:],'b-.')
        plt.plot(tau_range,Szt[2,:],'k-')
        plt.xlabel('Jt')
        plt.ylim(-0.6,0.6)
        plt.ylabel('$<S_{1}^{j}>/\hbar$')
        plt.legend(['j=x','j=y','j=z'])
        plt.show()


    def N(self, α, β, γ, circ, q1, q2):
        circ.rz(-0.5*np.pi,q2)
        circ.cx(q2,q1)
        circ.rz(0.5*np.pi-2*γ,q1)
        circ.ry(2.0*α-0.5*np.pi,q2)
        circ.cx(q1,q2)
        circ.ry(0.5*np.pi-2.0*β,q2)
        circ.cx(q2,q1)
        circ.rz(0.5*np.pi,q1)
    def clear_circuits(self):
        self.circuit_diagram_label_1.clear()
        self.circuit_diagram_label_2.clear()
        self.circuit_diagram_label_3.clear()
    def update_spin_system_parameters(self):
        self.selected_spin = self.spin_system_dropdown.currentIndex()
        params_to_show = {
            0: [0, 1, 2],  # Spin-1: Constant µ, Magnetic field strength, Time interval
            1: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11],  # Spin-2: Time interval, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2
            2: [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]  # Spin-3: Time interval, Jx, Jy, Jz, theta_1, phi_1, lambda_1, theta_2, phi_2, lambda_2, theta_3, phi_3, lambda_3, Number of Trotter steps
        }
        
        for i, (label, edit) in enumerate(zip(self.hamiltonian_labels, self.hamiltonian_edits)):
            if i in params_to_show[self.selected_spin]:
                label.show()
                edit.show()
            else:
                label.hide()
                edit.hide()        
    def onSubmit(self):
        # Gather input data
        circuit_params = {}
        for i, edit in enumerate(self.hamiltonian_edits):
            if edit.isVisible():
                circuit_params[self.hamiltonian_labels[i].text()] = edit.text()
        
        self.build_circuit(circuit_params, self.selected_spin)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = QuantumCircuitGUI()
    sys.exit(app.exec())
