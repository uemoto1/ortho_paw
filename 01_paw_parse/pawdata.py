#!/usr/bin/env python3
import gzip
import numpy as np
import xml.dom.minidom

class PAWData:

    def __init__(self, file_xml_gz):

        with gzip.open(file_xml_gz, "rt") as fh_gzip:
            xml_code = fh_gzip.read()

        root = xml.dom.minidom.parseString(xml_code)
        paw_dataset = root.getElementsByTagName("paw_dataset")[0]
        atom = paw_dataset.getElementsByTagName("atom")[0]

        self.atom_symbol = atom.getAttribute("symbol")
        self.atom_Z = float(atom.getAttribute("Z"))
        self.atom_core = float(atom.getAttribute("core"))
        self.atom_valence = float(atom.getAttribute("valence"))

        xc_functional = paw_dataset.getElementsByTagName("xc_functional")[0]
        self.xc_functional_type = xc_functional.getAttribute("type")
        self.xc_functional_name = xc_functional.getAttribute("name")

        generator = paw_dataset.getElementsByTagName("generator")[0]
        self.generator_type = generator.getAttribute("type")
        self.generator_name = generator.getAttribute("name")
         
        ae_energy = paw_dataset.getElementsByTagName("ae_energy")[0]
        self.ae_energy_kinetic = float(ae_energy.getAttribute("kinetic"))
        self.ae_energy_xc = float(ae_energy.getAttribute("xc"))
        self.ae_energy_electrostatic = float(ae_energy.getAttribute("electrostatic"))
        self.ae_energy_total = float(ae_energy.getAttribute("total"))

        core_energy = paw_dataset.getElementsByTagName("core_energy")[0]
        self.core_energy_kinetic = float(core_energy.getAttribute("kinetic"))

        paw_radius = paw_dataset.getElementsByTagName("paw_radius")[0]
        self.paw_radius_rc = float(paw_radius.getAttribute("rc"))
        
        valence_states = paw_dataset.getElementsByTagName("valence_states")[0]
        self.state_l = {}
        self.state_n = {}
        self.state_f = {}
        self.state_rc = {}
        self.state_e = {}
        for state in valence_states.getElementsByTagName("state"):
            state_id = state.getAttribute("id")
            self.state_l[state_id] = float(state.getAttribute("l"))
            self.state_e[state_id] = float(state.getAttribute("e"))
            self.state_rc[state_id] = float(state.getAttribute("rc"))
            state_n = state.getAttribute("n")
            if state_n.strip():
                self.state_n[state_id] = int(state_n)
            state_f = state.getAttribute("f")
            if state_f.strip():
                self.state_f[state_id] = float(state_f)

        radial_grid = paw_dataset.getElementsByTagName("radial_grid")[0]
        self.radial_grid_eq = radial_grid.getAttribute("eq")
        self.radial_grid_a = float(radial_grid.getAttribute("a"))
        self.radial_grid_d = float(radial_grid.getAttribute("d"))
        self.radial_grid_istart = int(radial_grid.getAttribute("istart"))
        self.radial_grid_iend = int(radial_grid.getAttribute("iend"))
        self.radial_grid_id = radial_grid.getAttribute("id")
        values = radial_grid.getElementsByTagName("values")[0]
        self.radial_grid_values = np.fromstring(values.childNodes[0].data, sep=" ")
        derivatives = radial_grid.getElementsByTagName("derivatives")[0]
        self.radial_grid_derivatives = np.fromstring(derivatives.childNodes[0].data, sep=" ")

        shape_function = paw_dataset.getElementsByTagName("shape_function")[0]
        self.shape_function_type = shape_function.getAttribute("type")
        self.shape_function_rc = float(shape_function.getAttribute("rc"))

        ae_core_density = paw_dataset.getElementsByTagName("ae_core_density")[0]
        self.ae_core_density_grid = ae_core_density.getAttribute("grid")
        self.ae_core_density_rc = float(ae_core_density.getAttribute("rc"))
        self.ae_core_density = np.fromstring(ae_core_density.childNodes[0].data, sep=" ")

        pseudo_core_density = paw_dataset.getElementsByTagName("pseudo_core_density")[0]
        self.pseudo_core_density_grid = pseudo_core_density.getAttribute("grid")
        self.pseudo_core_density_rc = float(pseudo_core_density.getAttribute("rc"))
        self.pseudo_core_density = np.fromstring(pseudo_core_density.childNodes[0].data, sep=" ")

        pseudo_valence_density = paw_dataset.getElementsByTagName("pseudo_valence_density")[0]
        self.pseudo_valence_density_grid = pseudo_valence_density.getAttribute("grid")
        self.pseudo_valence_density_rc = float(pseudo_valence_density.getAttribute("rc"))
        self.pseudo_valence_density = np.fromstring(pseudo_valence_density.childNodes[0].data, sep=" ")

        zero_potential = paw_dataset.getElementsByTagName("zero_potential")[0]
        self.zero_potential_grid = zero_potential.getAttribute("grid")
        self.zero_potential_rc = float(zero_potential.getAttribute("rc"))
        self.zero_potential = np.fromstring(zero_potential.childNodes[0].data, sep=" ")

        blochl_local_ionic_potential = paw_dataset.getElementsByTagName("blochl_local_ionic_potential")[0]
        self.blochl_local_ionic_potential_grid = blochl_local_ionic_potential.getAttribute("grid")
        self.blochl_local_ionic_potential_rc = float(blochl_local_ionic_potential.getAttribute("rc"))
        self.blochl_local_ionic_potential = np.fromstring(blochl_local_ionic_potential.childNodes[0].data, sep=" ")

        self.ae_partial_wave_grid = {}; self.ae_partial_wave = {}
        self.pseudo_partial_wave_grid = {}; self.pseudo_partial_wave = {}
        self.projector_function_grid = {}; self.projector_function = {}

        for ae_partial_wave in paw_dataset.getElementsByTagName("ae_partial_wave"):
            state = ae_partial_wave.getAttribute("state")
            self.ae_partial_wave_grid[state] = ae_partial_wave.getAttribute("grid")
            self.ae_partial_wave[state] = np.fromstring(ae_partial_wave.childNodes[0].data, sep=" ")

        for pseudo_partial_wave in paw_dataset.getElementsByTagName("pseudo_partial_wave"):
            state = pseudo_partial_wave.getAttribute("state")
            self.pseudo_partial_wave_grid[state] = pseudo_partial_wave.getAttribute("grid")
            self.pseudo_partial_wave[state] = np.fromstring(pseudo_partial_wave.childNodes[0].data, sep=" ")

        for projector_function in paw_dataset.getElementsByTagName("projector_function"):
            state = projector_function.getAttribute("state")
            self.projector_function_grid[state] = projector_function.getAttribute("grid")
            self.projector_function[state] = np.fromstring(projector_function.childNodes[0].data, sep=" ")

        kinetic_energy_differences = paw_dataset.getElementsByTagName("kinetic_energy_differences")[0]
        self.kinetic_energy_differences =  np.fromstring(kinetic_energy_differences.childNodes[0].data, sep=" ")

        exact_exchange_X_matrix = paw_dataset.getElementsByTagName("exact_exchange_X_matrix")[0]
        self.exact_exchange_X_matrix =  np.fromstring(exact_exchange_X_matrix.childNodes[0].data, sep=" ")


    def dump_vars(self):
        print("# atom_symbol = %s" % self.atom_symbol)
        print("# atom_Z = %s" % self.atom_Z)
        print("# atom_core = %s" % self.atom_core)
        print("# atom_valence = %s" % self.atom_valence)
        print("# xc_functional_type = %s" % self.xc_functional_type)
        print("# xc_functional_name = %s" % self.xc_functional_name)
        print("# generator_type = %s" % self.generator_type)
        print("# generator_name = %s" % self.generator_name)
        print("# ae_energy_kinetic = %s" % self.ae_energy_kinetic)
        print("# ae_energy_xc = %s" % self.ae_energy_xc)
        print("# ae_energy_electrostatic = %s" % self.ae_energy_electrostatic)
        print("# ae_energy_total = %s" % self.ae_energy_total)
        print("# core_energy_kinetic = %s" % self.core_energy_kinetic)
        print("# paw_radius_rc = %s" % self.paw_radius_rc)
        print("# state_l = %s" % self.state_l)
        print("# state_e = %s" % self.state_e)
        print("# state_rc = %s" % self.state_rc)
        print("# state_n = %s" % self.state_n)
        print("# state_f = %s" % self.state_f)
        print("# radial_grid_eq = %s" % self.radial_grid_eq)
        print("# radial_grid_a = %s" % self.radial_grid_a)
        print("# radial_grid_d = %s" % self.radial_grid_d)
        print("# radial_grid_istart = %s" % self.radial_grid_istart)
        print("# radial_grid_iend = %s" % self.radial_grid_iend)
        print("# radial_grid_id = %s" % self.radial_grid_id)
        print("# radial_grid_values = %s" % self.radial_grid_values)
        print("# radial_grid_derivatives = %s" % self.radial_grid_derivatives)
        print("# shape_function_type = %s" % self.shape_function_type)
        print("# shape_function_rc = %s" % self.shape_function_rc)
        print("# ae_core_density_grid = %s" % self.ae_core_density_grid)
        print("# ae_core_density_rc = %s" % self.ae_core_density_rc)
        print("# ae_core_density = %s" % self.ae_core_density)
        print("# pseudo_core_density_grid = %s" % self.pseudo_core_density_grid)
        print("# pseudo_core_density_rc = %s" % self.pseudo_core_density_rc)
        print("# pseudo_core_density = %s" % self.pseudo_core_density)
        print("# pseudo_valence_density_grid = %s" % self.pseudo_valence_density_grid)
        print("# pseudo_valence_density_rc = %s" % self.pseudo_valence_density_rc)
        print("# pseudo_valence_density = %s" % self.pseudo_valence_density)
        print("# zero_potential_grid = %s" % self.zero_potential_grid)
        print("# zero_potential_rc = %s" % self.zero_potential_rc)
        print("# zero_potential = %s" % self.zero_potential)
        print("# blochl_local_ionic_potential_grid = %s" % self.blochl_local_ionic_potential_grid)
        print("# blochl_local_ionic_potential_rc = %s" % self.blochl_local_ionic_potential_rc)
        print("# blochl_local_ionic_potential = %s" % self.blochl_local_ionic_potential)
        print("# ae_partial_wave_grid = %s" % self.ae_partial_wave_grid)
        print("# ae_partial_wave = %s" % self.ae_partial_wave)
        print("# pseudo_partial_wave_grid = %s" % self.pseudo_partial_wave_grid)
        print("# pseudo_partial_wave = %s" % self.pseudo_partial_wave)
        print("# projector_function_grid = %s" % self.projector_function_grid)
        print("# projector_function = %s" % self.projector_function)
        print("# kinetic_energy_differences = %s" % self.kinetic_energy_differences)
        print("# exact_exchange_X_matrix = %s" % self.exact_exchange_X_matrix)
        print("# zero_potential = %s" % self.zero_potential)





