#PBUTs SIMULATOR - MODELLO FORWARD EULER (HDF POST-DILUIZIONE) + LETTURA AUTOMATICA DA EXCEL
#Script integrato: unisce ForwardEuler_Maheshwari.py e Automatizz_excel.py
#06DICEMBRE2025 - Integrazione

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime
import threading
import os
import pandas as pd

# =====================================================================
# 1. PARAMETRI DI DEFAULT DEL MODELLO
# =====================================================================

Filter_models = {
    "Fx100": {
        "KoA": 4.08305e-4,  # m^3/min valore preso dai dati KoA per red phenol senza BSA e diminuito del 20%
        "Lfiber": 22.8e-2,     # m
        "N": 16896,     # fibers dato preso dal dropbox
        "A": 2.21e-8,   # m^2 fiber area
        "Ad": 4.30e-8,  # m^2 fiber area (dialysate side)
        "S_T": 1
    },

    "Solacea21H": {
        "KoA": 5.42956e-4,  # m^3/min valore preso dai dati KoA per red phenol senza BSA e diminuito del 20%
        "Lfiber": 25.4e-2,     # m
        "N": 12960,     # fibers
        "A": 2.78e-8,   # m^2 fiber area
        "Ad": 4.45e-8,   # m^2 fiber area (dialysate side)
        "S_T": 1
    }
}

Toxin_parameters = {
    # indoxyl sulfate IS
    "IS": {
        "MwT": 213.21,  # g/mol
        "k1": 3.6e4,  # m^3/(mol·min)  controllato! valore preso da Protein-bound uremic toxin kinetics modeling for a new dialysis device
        "G": 0.0285 / 213.21 * 1e-3,  # mol/min         0.0285 mg/min
        "Kip_T": 1210e-6,  # m^3/min
        "Kic_T": 103e-6,   # m^3/min
    },
    # P-cresol sulfate pCS
    "pCS": {
        "MwT": 188.22,  # g/mol
        "k1": 1e4,  # m^3/(mol·min)
        "G": 0.0277 / 188.22 * 1e-3,  # mol/min
        "Kip_T": 1060e-6,  # m^3/min
        "Kic_T": 98e-6   # m^3/min
    }
}

Model_Parameters = {
    # --- HDF post venous tube parameters ---#
    "n_segments_dialyzer": 40,
    "n_segments_tube": 100,
    "dt": 0.00001,  # Time step for the simulation (min)
    "total_time": 240,  # Total time for the simulation (min)
    "out_time": 1,  # Time interval for output (min)
    "MwP": 66.1 * 1e3  # Molecular weight of albumin protein (g/mol)
}

Dialyzation_settings = {
    "Qp": 294e-6,      # Plasma flow (m^3/min).
    "Quf": 95.25862e-6,  # Ultrafiltration rate (m^3/min)
    "Qd": 417.67e-6,      # Dialysate flow (m^3/min)
    "Qr": 82.3276e-6       # Replacement fluid rate (m^3/min)
}

Patient = {
    "Vpl0": 3.775e-3,  # Plasma volume in m^3
    "Vis0": 11.325e-3,  # Interstitial volume in m^3
    "Vic": 18.4e-3,  # Intracellular volume in m^3
    "Tpl0": 2.357 / 213.21,  # Initial plasma concentration in mol/m^3.
    "PTpl0": (31.02 - 2.357) / (213.21),  # Initial plasma toxin concentration in mol/m^3
    "Ppl0": 0.6,  # Initial plasma protein concentration in mol/m^3.  4g/dL -> 0.6 mol/m^3
    "Tis0": 2.357 / 213.21,  # Initial interstitial concentration in mol/m^3
    "PTis0": (31.02 - 2.357) / (2 * 213.21),  # Initial interstitial toxin concentration in mol/m^3     =la metà di PTpl0
    "Pis0": 0.3,  # Initial interstitial protein concentration in mol/m^3       =la metà di Ppl0
    "Tic0": 2.357 / 213.21  # Initial intracellular concentration in mol/m^3            T uguale in tutti i compartimenti
}

Tube_Parameters = {
    "Ltube": 0.5,  # Length of the HDF post venous tube (m)
    "Atube": 1e-3 ** 2 * 3.14,  # Cross-sectional area of the tube (m^2)
}


# =====================================================================
# 2. CLASSE DEL MODELLO (Forward Euler)
# =====================================================================

class PbutsPredictor:
    def __init__(self, filter_model="Fx100", Toxin="IS", Patient=Patient, Model_Parameters=Model_Parameters,
                 Tube_Parameters=Tube_Parameters, Dialyzation_settings=Dialyzation_settings, out_name="",
                 out_path="."):
        # load the parameters from the Toxin database
        # Parameters for the patient
        # Mw, k1, G, Kip_T, Kic_T
        for key, value in Filter_models[filter_model].items():
            setattr(self, key, value)
        for key, value in Toxin_parameters[Toxin].items():
            setattr(self, key, value)
        # Parameters for the patient
        for key, value in Patient.items():
            setattr(self, key, value)
        # Model parameters
        for key, value in Model_Parameters.items():
            setattr(self, key, value)
        # Dialyzer parameters
        for key, value in Dialyzation_settings.items():
            setattr(self, key, value)
        # Patient parameters
        for key, value in Patient.items():
            setattr(self, key, value)
        # Tube parameters
        for key, value in Tube_Parameters.items():
            setattr(self, key, value)

        self.out_path = out_path
        self.out_name = out_name
        self.MwPT = self.MwP + self.MwT  # Molecular weight of the complex (g/mol)
        self.dx = self.Lfiber / self.n_segments_dialyzer  # Spatial step size for the dialyzer
        self.total_time_steps = int(self.total_time / self.dt)
        self.out_steps = int(self.out_time / self.dt)
        self.dx_tube = self.Ltube / self.n_segments_tube  # Spatial step size for the tube
        self.out_steps = int(self.out_time / self.dt)  # Number of time steps to save results
        self.Qtube = self.Qp - self.Quf + self.Qr  # Total flow rate in the tube (m^3/min)
        # Check stability condition
        if self.dt / (self.dx ** 2) > 1:
            print(f'Warning: The stability condition is not satisfied (dt/(dx^2)={self.dt / (self.dx ** 2):.2f} > 1). The simulation may be unstable.')
        else:
            print(f'Stability condition is satisfied (dt/(dx^2)={self.dt / (self.dx ** 2):.2f} <= 1). The simulation should be stable.')
        # --- Initialize arrays for concentrations ---#
        if self.dt / (self.dx_tube ** 2) > 1:
            print(f'Warning: The stability condition is not satisfied (dt/(dx^2)={self.dt / (self.dx_tube ** 2):.2f} > 1). The simulation may be unstable.')
        else:
            print(f'Stability condition is satisfied (dt/(dx^2)={self.dt / (self.dx_tube ** 2):.2f} <= 1). The simulation should be stable.')

        # Peclet number and flux adjustment
        self.Pe = (self.Quf) / self.KoA
        self.Qdx = lambda x: self.Qd + (self.Lfiber - x) / self.Lfiber * (self.Quf)  # Dialyzer flow rate as a function of position
        self.Qpx = lambda x: self.Qp - (x) / self.Lfiber * (self.Quf)
        self.alpha_t = lambda Vis, Vpl: (Vis) / (Vis + Vpl)
        # print all the parameters
        print("Dialyzer model parameters:")
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")
        self.running = False

    def start(self):
        # Start the simulation in a separate thread
        if not self.running:
            self.running = True
            threading.Thread(target=self.run, daemon=True).start()

    def stop(self):
        self.running = False

    def run(self):
        # Initialize arrays for results
        Vpl = np.zeros(self.total_time_steps)
        Vis = np.zeros(self.total_time_steps)
        T_pl = np.zeros(self.total_time_steps)
        PT_pl = np.zeros(self.total_time_steps)
        P_pl = np.zeros(self.total_time_steps)
        T_is = np.zeros(self.total_time_steps)
        PT_is = np.zeros(self.total_time_steps)
        P_is = np.zeros(self.total_time_steps)
        T_ic = np.zeros(self.total_time_steps)
        T_dpl = np.zeros((2, self.n_segments_dialyzer))
        PT_dpl = np.zeros((2, self.n_segments_dialyzer))
        P_dpl = np.zeros((2, self.n_segments_dialyzer))
        T_d = np.zeros((2, self.n_segments_dialyzer))

        T_tpl = np.zeros((2, self.n_segments_tube))  # Plasma concentration in venous tube
        PT_tpl = np.zeros((2, self.n_segments_tube))  # Plasma toxin concentration in venous tube
        P_tpl = np.zeros((2, self.n_segments_tube))  # Plasma protein concentration in venous tube

        # Set initial conditions
        Vpl[0] = self.Vpl0
        Vis[0] = self.Vis0
        Vic = self.Vic  # Assuming constant intracellular volume
        T_pl[0] = self.Tpl0  # Convert to mol/m^3
        PT_pl[0] = self.PTpl0
        P_pl[0] = self.Ppl0
        T_is[0] = self.Tis0
        PT_is[0] = self.PTis0
        P_is[0] = self.Pis0
        T_ic[0] = self.Tic0
        T_dpl[0, :] = 0
        PT_dpl[0, :] = 0
        P_dpl[0, :] = 0
        T_d[0, :] = 0
        self.k2 = self.k1 * (self.Ppl0 * self.Tpl0) / (self.PTpl0)  # Dissociation rate constant (min^-1)

        j_dialyzer = np.arange(1, self.n_segments_dialyzer)  # Indices for the interior points
        k_dialyzer = j_dialyzer - 1  # Index for the dialyzer segment

        j_tube = np.arange(1, self.n_segments_tube)  # Indices for the interior points

        start_time = datetime.today().strftime('%Y%m%d_%H%M')

        self.out_results_patient = os.path.join(self.out_path, f'FE_results_patient_{self.out_name}_IS_{start_time}.txt')
        self.out_results_dialyzer = os.path.join(self.out_path, f'FE_results_dialyzer_{self.out_name}_IS_{start_time}.txt')
        self.out_results_dialyzer_kinetics = os.path.join(self.out_path, f'FE_dialyzer_kinetics_{self.out_name}_IS_{start_time}.txt')
        self.out_results_tube_kinetics = os.path.join(self.out_path, f'FE_tube_kinetics_{self.out_name}_IS_{start_time}.txt')

        with open(self.out_results_dialyzer, 'w') as dialyzer_file, \
             open(self.out_results_dialyzer_kinetics, 'w') as dialyzer_kinetics_file, \
             open(self.out_results_tube_kinetics, 'w') as tube_kinetics_file, \
             open(self.out_results_patient, 'w') as patient_file:
            patient_file.write('Time (min), T_pl (mg/L), PT_pl (mg/L), P_pl (mg/L), T_is (mg/L), PT_is (mg/L), P_is (mg/L), T_ic (mg/L)\n')
            dialyzer_file.write('Time (min), T_dpl_out (mg/L), PT_dpl_out (mg/L), P_dpl_out (mg/L), T_d_out (mg/L)\n')
            dialyzer_kinetics_file.write('Time (min), T_dpl_out (mg/L), PT_dpl_out (mg/L), P_dpl_out (mg/L), T_d_out (mg/L)\n')
            tube_kinetics_file.write('Time (min), T_tpl_out (mg/L), PT_tpl_out (mg/L), P_tpl_out (mg/L)\n')
            for i in tqdm(range(0, self.total_time_steps - 1)):
                if np.isnan(T_pl[i]):
                    print(f'Warning: NaN detected at time step {i}. exit.')
                    break
                elif not self.running:
                    print("Simulation stopped by user.")
                    break
                if i % self.out_steps == 0:
                    # Save results to files every out_time minutes
                    time_min = i * self.dt
                    dialyzer_file.write(f'{time_min:.2f}, {T_dpl[0, -1] * self.MwT:.5e}, {PT_dpl[0, -1] * self.MwT:.5e}, {P_dpl[0, -1] * self.MwP:.5e}, {T_d[0, 0] * self.MwT:.5e}\n')
                    dialyzer_file.flush()
                    dialyzer_kinetics_file.write(f'{time_min:.2f}, {(T_dpl[0, :] * self.MwT).tolist()}, {(PT_dpl[0, :] * self.MwT).tolist()}, {(P_dpl[0, :] * self.MwP).tolist()}, {(T_d[0, :] * self.MwT).tolist()}\n')
                    dialyzer_kinetics_file.flush()
                    tube_kinetics_file.write(f'{time_min:.2f}, {(T_tpl[0, :] * self.MwT).tolist()}, {(PT_tpl[0, :] * self.MwT).tolist()}, {(P_tpl[0, :] * self.MwP).tolist()}\n')
                    tube_kinetics_file.flush()
                    patient_file.write(f'{time_min:.2f}, {T_pl[i] * self.MwT:.5e}, {PT_pl[i] * self.MwT:.5e}, {P_pl[i] * self.MwP:.5e}, {T_is[i] * self.MwT:.5e}, {PT_is[i] * self.MwT:.5e}, {P_is[i] * self.MwP:.5e}, {T_ic[i] * self.MwT:.5e}\n')
                    patient_file.flush()

                # Calculate the alpha value based on the current volumes
                alpha = self.alpha_t(Vis[i], Vpl[i])  # Calculate alpha based on current volumes
                # Update volumes
                dVpl = -(1 - alpha) * self.Quf + self.Qr
                dVis = -alpha * self.Quf
                Vpl[i + 1] = Vpl[i] + dVpl * self.dt
                Vis[i + 1] = Vis[i] + dVis * self.dt
                # Assuming constant volume for intracellular water
                # Patient kinetics
                # Update the Toxin concentrations (T), Proteins (P), and their complexes (PT) in the plasma
                dT_pl = (-self.Qp * T_pl[i] + (self.Qp - self.Quf + self.Qr) * T_tpl[0, -1] + self.Kip_T * (T_is[i] - T_pl[i]) + alpha * self.Quf * T_is[i] + (self.k2 * PT_pl[i] - self.k1 * P_pl[i] * T_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * T_pl[i] / Vpl[i]
                T_pl[i + 1] = T_pl[i] + self.dt * dT_pl
                dPT_pl = (-self.Qp * PT_pl[i] + (self.Qp - self.Quf + self.Qr) * PT_tpl[0, -1] + (self.k1 * P_pl[i] * T_pl[i] - self.k2 * PT_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * PT_pl[i] / Vpl[i]
                PT_pl[i + 1] = PT_pl[i] + self.dt * dPT_pl
                dP_pl = (-self.Qp * P_pl[i] + (self.Qp - self.Quf + self.Qr) * P_tpl[0, -1] + (-self.k1 * P_pl[i] * T_pl[i] + self.k2 * PT_pl[i]) * Vpl[i]) / Vpl[i] - dVpl * P_pl[i] / Vpl[i]
                P_pl[i + 1] = P_pl[i] + self.dt * dP_pl
                # Update the Toxin concentrations (T), Proteins (P), and their complexes (PT) in the interstitial space
                dT_is = (self.Kic_T * (T_ic[i] - T_is[i]) - self.Kip_T * (T_is[i] - T_pl[i]) - alpha * self.Quf * T_is[i] + (self.k2 * PT_is[i] - self.k1 * P_is[i] * T_is[i]) * Vis[i]) / Vis[i] - dVis * T_is[i] / Vis[i]
                T_is[i + 1] = T_is[i] + self.dt * dT_is
                dPT_is = (self.k1 * P_is[i] * T_is[i] - self.k2 * PT_is[i]) - dVis * PT_is[i] / Vis[i]
                PT_is[i + 1] = PT_is[i] + self.dt * dPT_is
                dP_is = (-self.k1 * P_is[i] * T_is[i] + self.k2 * PT_is[i]) - dVis * P_is[i] / Vis[i]
                P_is[i + 1] = P_is[i] + self.dt * dP_is
                # Update the Toxin concentrations (T) in the intracellular space, assuming that the Proteins (P) and thus their complexes (PT) is equal to 0
                dT_ic = (self.G - self.Kic_T * (T_ic[i] - T_is[i])) / Vic  # intracellular volume constant
                T_ic[i + 1] = T_ic[i] + self.dt * dT_ic

                # Boundary conditions
                T_dpl[1, 0] = T_pl[i + 1]  # Initial condition for the dialyzer Toxin concentration
                PT_dpl[1, 0] = PT_pl[i + 1]  # Initial condition for the dialyzer complex concentration
                P_dpl[1, 0] = P_pl[i + 1]  # Initial condition for the dialyzer Protein concentration
                T_d[1, -1] = 0  # Initial condition for the dialysate Toxin concentration

                # Dialyzer kinetics
                # Update the Toxin concentrations (T), Proteins (P), and their complexes (PT) in the dialyzer
                dx_Tpl = (T_dpl[0, j_dialyzer] - T_dpl[0, j_dialyzer - 1]) / self.dx  # First order finite difference for T_dpl
                dT_dpl = (-1 / (self.N * self.A) * (dx_Tpl * self.Qpx(j_dialyzer * self.dx) - (self.Quf) / self.Lfiber * T_dpl[0, j_dialyzer]) - 1 / (self.N * self.A * self.Lfiber) * (self.KoA * (T_dpl[0, j_dialyzer] - T_d[0, j_dialyzer])) * self.Pe / (np.exp(self.Pe) - 1) + 1 / (self.N * self.A * self.Lfiber) * (-((self.Quf) * T_dpl[0, j_dialyzer] * self.S_T)) + (-self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] + self.k2 * PT_dpl[0, j_dialyzer]))
                T_dpl[1, j_dialyzer] = T_dpl[0, j_dialyzer] + self.dt * dT_dpl

                dx_PTdpl = (PT_dpl[0, j_dialyzer] - PT_dpl[0, j_dialyzer - 1]) / self.dx  # First order finite difference for PT_dpl
                dPT_dpl = (-1 / (self.N * self.A) * (dx_PTdpl * self.Qpx(j_dialyzer * self.dx) - (self.Quf) / self.Lfiber * PT_dpl[0, j_dialyzer]) + (self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] - self.k2 * PT_dpl[0, j_dialyzer]))
                PT_dpl[1, j_dialyzer] = PT_dpl[0, j_dialyzer] + self.dt * dPT_dpl

                dx_Pdpl = (P_dpl[0, j_dialyzer] - P_dpl[0, j_dialyzer - 1]) / self.dx  # First order finite difference for P_dpl
                dP_dpl = (-1 / (self.N * self.A) * (dx_Pdpl * self.Qpx(j_dialyzer * self.dx) - (self.Quf) / self.Lfiber * P_dpl[0, j_dialyzer]) + (-self.k1 * P_dpl[0, j_dialyzer] * T_dpl[0, j_dialyzer] + self.k2 * PT_dpl[0, j_dialyzer]))
                P_dpl[1, j_dialyzer] = P_dpl[0, j_dialyzer] + self.dt * dP_dpl

                # Update the Toxin concentrations (T) in the dialysate
                dx_Td = (T_d[0, k_dialyzer + 1] - T_d[0, k_dialyzer - 1]) / (2 * self.dx)  # First order finite difference for T_d
                dx_Td[0] = (T_d[0, 1] - T_d[0, 0]) / (self.dx)  # First order finite difference for T_pl
                dT_d = 1 / (self.N * self.Ad) * (dx_Td * self.Qdx(self.dx * k_dialyzer) - (self.Quf) / self.Lfiber * T_d[0, k_dialyzer]) + 1 / (self.N * self.Ad * self.Lfiber) * (self.KoA * (T_dpl[0, k_dialyzer] - T_d[0, k_dialyzer])) * self.Pe / (np.exp(self.Pe) - 1) - 1 / (self.N * self.Ad * self.Lfiber) * (-((self.Quf) * T_dpl[0, k_dialyzer] * self.S_T))
                T_d[1, k_dialyzer] = T_d[0, k_dialyzer] + self.dt * dT_d

                T_dpl[0, :] = T_dpl[1, :]  # Update the dialysate concentration for the next iteration
                PT_dpl[0, :] = PT_dpl[1, :]  # Update the dialysate complex concentration for the next iteration
                P_dpl[0, :] = P_dpl[1, :]  # Update the dialysate protein concentration for the next iteration
                T_d[0, :] = T_d[1, :]  # Update the dialysate concentration for the next iteration

                # To add postra hemodiafiltration
                # Boundary conditions for the postra hemodiafiltration

                T_tpl[1, 0] = T_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube
                PT_tpl[1, 0] = PT_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube
                P_tpl[1, 0] = P_dpl[1, -1] * self.Qpx(self.Lfiber) / self.Qtube  # Initial condition for the plasma protein concentration in the venous tube

                dx_tTpl = (T_tpl[0, j_tube] - T_tpl[0, j_tube - 1]) / (self.dx_tube)  # First order finite difference for T_tpl
                dt_tTpl = -(self.Qtube / self.Atube) * dx_tTpl + (-self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] + self.k2 * PT_tpl[0, j_tube])
                T_tpl[1, j_tube] = T_tpl[0, j_tube] + self.dt * dt_tTpl

                dx_tPTpl = (PT_tpl[0, j_tube] - PT_tpl[0, j_tube - 1]) / (self.dx_tube)  # First order finite difference for PT_tpl
                dt_tPTpl = -(self.Qtube / self.Atube) * dx_tPTpl + (self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] - self.k2 * PT_tpl[0, j_tube])
                PT_tpl[1, j_tube] = PT_tpl[0, j_tube] + self.dt * dt_tPTpl

                dx_tPpl = (P_tpl[0, j_tube] - P_tpl[0, j_tube - 1]) / (self.dx_tube)  # First order finite difference for P_ppl
                dt_tPpl = -(self.Qtube / self.Atube) * dx_tPpl + (-self.k1 * P_tpl[0, j_tube] * T_tpl[0, j_tube] + self.k2 * PT_tpl[0, j_tube])
                P_tpl[1, j_tube] = P_tpl[0, j_tube] + self.dt * dt_tPpl

                T_tpl[0, :] = T_tpl[1, :]  # Update the plasma concentration in the venous tube for the next iteration
                PT_tpl[0, :] = PT_tpl[1, :]  # Update the plasma toxin concentration in the venous tube for the next iteration
                P_tpl[0, :] = P_tpl[1, :]  # Update the plasma protein concentration in the venous tube for the next iteration
                # Print the results for debugging
        self.running = False


# =====================================================================
# 3. FUNZIONE DI AUTOMAZIONE: LEGGE IL DATABASE EXCEL E LANCIA LE SIMULAZIONI
# =====================================================================

def avvia_simulazione_mirata(
    excel_path="database_patients_IS_T6.xlsx",
    pazienti_target=None,
    out_path="./results_IS_T6",
    toxin="IS",
):
    """
    Legge il database Excel dei pazienti e lancia una simulazione PbutsPredictor
    per ciascun paziente indicato in pazienti_target (o per tutti se None).

    Parametri
    ---------
    excel_path : str
        Percorso del file Excel con i dati dei pazienti (una riga per paziente).
    pazienti_target : list[str] o None
        Lista degli ID_Paziente da simulare. Se None, vengono simulati tutti
        i pazienti presenti nel file.
    out_path : str
        Cartella in cui salvare i risultati (creata se non esiste).
    toxin : str
        Tossina da simulare ("IS" o "pCS"), deve corrispondere alle colonne
        del file Excel.
    """
    if pazienti_target is None:
        pazienti_target = ["M01"]

    os.makedirs(out_path, exist_ok=True)

    df = pd.read_excel(excel_path)

    # Filtro il dataframe: prendo solo le righe dove ID_Paziente è nella lista
    df_filtrato = df[df['ID_Paziente'].isin(pazienti_target)]

    # Peso molecolare per la conversione
    mw_toxin = Toxin_parameters[toxin]["MwT"]

    for index, row in df_filtrato.iterrows():
        print(f"Esecuzione per ID_Paziente = {row['ID_Paziente']}")

        # Preparazione dizionari (stessa logica dello script originale)
        dati_paziente = {
            "Vpl0": row['Vpl0'] * 1e-3,
            "Vis0": row['Vis0'] * 1e-3,
            "Vic": row['Vic'] * 1e-3,
            "Tpl0": row['Tpl0'] / mw_toxin,
            "PTpl0": row['PTpl0'] / mw_toxin,
            "Ppl0": row['Ppl0'] * 0.15045,
            "Tis0": row['Tis0'] / mw_toxin,
            "PTis0": row['PTis0'] / mw_toxin,
            "Pis0": row['Pis0'] * 0.15045,
            "Tic0": row['Tic0'] / mw_toxin
        }

        settaggi_dialisi = {
            "Qp": row['Qp'] * 1e-6,
            "Quf": row['Quf'] * 1e-6,
            "Qd": row['Qd'] * 1e-6,
            "Qr": row['Qr'] * 1e-6,
        }

        # Istanza del simulatore
        print(f"Filtro utilizzato: {row['Modello_Filtro']}")
        sim = PbutsPredictor(
            filter_model=row['Modello_Filtro'],
            Toxin=toxin,
            Patient=dati_paziente,
            Dialyzation_settings=settaggi_dialisi,
            out_name=row['ID_Paziente'],
            out_path=out_path
        )

        # Esecuzione
        print(f"Inizio simulazione per {row['ID_Paziente']}...")
        sim.running = True
        sim.run()
        print(f"Simulazione per {row['ID_Paziente']} completata!")


if __name__ == "__main__":
    # --- MODIFICA QUI I PARAMETRI DELLA CAMPAGNA DI SIMULAZIONE ---
    avvia_simulazione_mirata(
        excel_path="database_patients_IS_T6.xlsx",
        pazienti_target=["M01"],   # es. ["M01", "M02", ...] oppure None per tutti
        out_path="./results_IS_T6",
        toxin="IS",
    )
