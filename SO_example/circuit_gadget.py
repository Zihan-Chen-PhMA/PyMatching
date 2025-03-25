from __future__ import annotations
from surface_codes import *
from typing import Callable
import copy
import pymatching
      


class Meta_Circuit():

    def __init__(self, code_list: list[Union[Rotated_surface_code]], 
                 circuit_filename: str,
                 noise_probability: float = 0, 
                 control_functional: Callable = None) -> None:
        self.code_list = code_list
        self.code_name_dict: dict[str,Union[Rotated_surface_code]] = {}
        for code in self.code_list:
            self.code_name_dict[code.code_name] = code
        self.qubit_collection = []
        self.qubit_network: dict[tuple,Qubit] = {}
        self._collect_qubits()
        self.code_name_SE_instruction: dict[str,SE_instruction_set] = {}
        for code_name in self.code_name_dict:
            self.code_name_SE_instruction[code_name] = self._SE_compile(code_name)

        self.circuit_filename = circuit_filename
        self.circ = Circuit_helper(self.qubit_collection,
                            circuit_filename)
        self.noise_probability = noise_probability
        self.growth_time = -1
        self.circ.initialize()
        self.logic_qubits = []

    def _collect_qubits(self) -> None:
        self.qubit_collection = []
        id_temp = 0
        for code in self.code_list:
            for qubit in code.qubit_network.values():
                if qubit.pos not in self.qubit_network:
                    qubit_new = Qubit(id_temp,qubit.pos)
                    self.qubit_collection.append(qubit_new)
                    self.qubit_network[qubit.pos] = qubit_new

    def _SE_compile(self, code_name: str) -> SE_instruction_set:
        code_temp = self.code_name_dict[code_name]
        z_ctrl_dict: dict[int,list] = {}
        z_targ_dict: dict[int,list] = {}
        x_ctrl_dict: dict[int,list] = {}
        x_targ_dict: dict[int,list] = {}
        mid_cycle_steps = []
        
        for qubit in code_temp.z_check_collection + code_temp.x_check_collection:
            for step_str in qubit.neighbor:
                if int(step_str) not in mid_cycle_steps:
                    mid_cycle_steps.append(int(step_str))
        mid_cycle_steps = sorted(mid_cycle_steps)

        for step in mid_cycle_steps:
            z_ctrl_dict[step] = []
            z_targ_dict[step] = []
            x_ctrl_dict[step] = []
            x_targ_dict[step] = []
            for qubit in code_temp.z_check_collection:
                if str(step) in qubit.neighbor:
                    ctrl_pos = qubit.neighbor[str(step)].pos
                    targ_pos = qubit.pos
                    z_ctrl_dict[step].append(self.qubit_network[ctrl_pos])
                    z_targ_dict[step].append(self.qubit_network[targ_pos])
            for qubit in code_temp.x_check_collection:
                if str(step) in qubit.neighbor:
                    ctrl_pos = qubit.pos
                    targ_pos = qubit.neighbor[str(step)].pos
                    x_ctrl_dict[step].append(self.qubit_network[ctrl_pos])
                    x_targ_dict[step].append(self.qubit_network[targ_pos])
        
        ret_se_instruction = SE_instruction_set(code_name, mid_cycle_steps,
                                                z_ctrl_dict, z_targ_dict,
                                                x_ctrl_dict, x_targ_dict)
        return ret_se_instruction


    def SE_round(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                        noise=single_qubit_noise)
        self.circ.single_qubit_gate('I', self.qubit_collection,
                                    noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)
        for check in x_checks_in_circ + z_checks_in_circ:
            self.circ.detector([check], self.circ.time_temp - 1,
                               [check], self.circ.time_temp,
                        (check.pos[0],check.pos[1], self.circ.time_temp))

        self.circ.clock_plus_one()
        self.circ.tick()
        
        pass

    def SE_round_circ(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        x_checks_in_circ: list[Qubit] = []
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        z_checks_in_circ: list[Qubit] = []
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])

        self.circ.reset(x_checks_in_circ + z_checks_in_circ,
                        noise = MR_noise)
        self.circ.single_qubit_gate('H', x_checks_in_circ)
        self.circ.single_qubit_gate('I', self.qubit_collection,
                                    noise = single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)


    def X_init(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
            R_noise = ['X_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
            R_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        data_qubits_in_circ: list[Qubit] = []
        z_checks_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        
        self.circ.reset(data_qubits_in_circ + z_checks_in_circ \
                        + x_checks_in_circ, 
                        noise = R_noise)
        self.circ.single_qubit_gate('H', x_checks_in_circ + data_qubits_in_circ,
                                         noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)
        for qubit in x_checks_in_circ:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass

    def Z_init(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            two_qubit_noise = ['DEPOLARIZE2',self.noise_probability]
            MR_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
            R_noise = ['X_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            two_qubit_noise = None
            MR_noise = None
            MX_noise = None
            R_noise = None
        code_temp = self.code_name_dict[code_name]
        se_instruction = self.code_name_SE_instruction[code_name]
        data_qubits_in_circ: list[Qubit] = []
        z_checks_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for z_check in code_temp.z_check_collection:
            z_checks_in_circ.append(self.qubit_network[z_check.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])
        
        self.circ.reset(data_qubits_in_circ + z_checks_in_circ \
                        + x_checks_in_circ, 
                        noise = MR_noise)
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                         noise=single_qubit_noise)
        self.circ.tick()
        for i in se_instruction.mid_cycle_steps:
            ctrl_list = []
            targ_list = []
            try:
                ctrl_list += se_instruction.z_ctrl_dict[i]
                targ_list += se_instruction.z_targ_dict[i]
            except:
                pass
            try:
                ctrl_list += se_instruction.x_ctrl_dict[i]
                targ_list += se_instruction.x_targ_dict[i]
            except:
                pass

            self.circ.two_qubit_gate('CX', ctrl_list,
                                           targ_list,
                                           noise=two_qubit_noise)
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', x_checks_in_circ,
                                    noise=single_qubit_noise)
        self.circ.tick()
        self.circ.measure_reset(x_checks_in_circ + z_checks_in_circ, 
                                noise = MR_noise)
        for qubit in z_checks_in_circ:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass

    def X_meas(self, code_name: str, cleaness: bool = False) -> None:
        if cleaness == False:
            single_qubit_noise = ['DEPOLARIZE1',self.noise_probability]
            MZ_noise = ['X_ERROR', self.noise_probability]
            MX_noise = ['Z_ERROR', self.noise_probability]
        else:
            single_qubit_noise = None
            MZ_noise = None
            MX_noise = None
        code_temp = self.code_name_dict[code_name]
        data_qubits_in_circ: list[Qubit] = []
        x_checks_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])
        for x_check in code_temp.x_check_collection:
            x_checks_in_circ.append(self.qubit_network[x_check.pos])

        self.circ.single_qubit_gate('H',data_qubits_in_circ,
                                    noise= single_qubit_noise)
        self.circ.measure(data_qubits_in_circ,
                          noise = MZ_noise)
        for x_check in code_temp.x_check_collection:
            check_pos = x_check.pos
            x_check_in_circ = self.qubit_network[check_pos]
            neighbors_in_circ = []
            for neighbor in x_check.neighbor.values():
                neighbors_in_circ.append(self.qubit_network[neighbor.pos])
            self.circ.detector([x_check_in_circ], self.circ.time_temp - 1,
                        neighbors_in_circ,
                        self.circ.time_temp,
                        (x_check_in_circ.pos[0],x_check_in_circ.pos[1], self.circ.time_temp))
        
        logic_x_in_circ = []
        for qubit in code_temp.logic_x_collection:
            logic_x_in_circ.append(self.qubit_network[qubit.pos])
        self.circ.observable(logic_x_in_circ,self.circ.time_temp)

        pass

    def Z_meas(self, code_name: str) -> None:
        code_temp = self.code_name_dict[code_name]
        data_qubits_in_circ: list[Qubit] = []
        for data in code_temp.data_qubits_collection:
            data_qubits_in_circ.append(self.qubit_network[
                                            data.pos])

        self.circ.measure(data_qubits_in_circ,
                          noise = ['X_ERROR', self.noise_probability])
        for z_check in code_temp.z_check_collection:
            check_pos = z_check.pos
            z_check_in_circ = self.qubit_network[check_pos]
            neighbors_in_circ = []
            for neighbor in z_check.neighbor.values():
                neighbors_in_circ.append(self.qubit_network[neighbor.pos])
            self.circ.detector([z_check_in_circ], self.circ.time_temp - 1,
                        neighbors_in_circ,
                        self.circ.time_temp,
                        (z_check_in_circ.pos[0],z_check_in_circ.pos[1], self.circ.time_temp))
        
        logic_z_in_circ = []
        for qubit in code_temp.logic_z_collection:
            logic_z_in_circ.append(self.qubit_network[qubit.pos])
        self.circ.observable(logic_z_in_circ,self.circ.time_temp)
        pass


    def check_circuit_distance(self) -> int:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        explained_error = circuit.search_for_undetectable_logical_errors(
                            dont_explore_detection_event_sets_with_size_above=6,
                            dont_explore_edges_with_degree_above=3,
                            dont_explore_edges_increasing_symptom_degree=False)
        print(len(explained_error))
        for item in explained_error:
            print(item)
        return len(explained_error)
    

        
    def get_stim_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        return circuit
   

class SE_instruction_set():

    def __init__(self, code_name: str,
                 mid_cycle_steps: list,
                 z_ctrl_dict: dict[int,list[Qubit]],
                 z_targ_dict: dict[int,list[Qubit]],
                 x_ctrl_dict: dict[int,list[Qubit]],
                 x_targ_dict: dict[int,list[Qubit]]
                 ):
        
        self.code_name = code_name
        self.mid_cycle_steps = mid_cycle_steps
        self.z_ctrl_dict = z_ctrl_dict
        self.z_targ_dict = z_targ_dict
        self.x_ctrl_dict = x_ctrl_dict
        self.x_targ_dict = x_targ_dict



