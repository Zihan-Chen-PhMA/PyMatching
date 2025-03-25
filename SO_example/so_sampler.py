from __future__ import annotations
import collections
import dataclasses
import heapq
import math
import time
from typing import Literal, cast, Any, AbstractSet, Union, Optional

import numpy as np
import pymatching
import sinter
import stim
from dem_parsor import *


class SOXSampler(sinter.Sampler):
    def compiled_sampler_for_task(self, task: sinter.Task) -> CompiledSOXSampler:
        return CompiledSOXSampler.from_task(task)
    pass

class CompiledSOXSampler(sinter.CompiledSampler):
    def __init__(self, dem: stim.DetectorErrorModel,
                 circuit: stim.Circuit,
                 post_select_det_ids: frozenset,
                 decoder: Optional[pymatching.Matching] = None):
        self.dem = dem
        self.num_detectors = dem.num_detectors
        self.circuit = circuit
        self.sampler = self.circuit.compile_detector_sampler()
        self.post_select_det_ids = post_select_det_ids
        self.decoder = decoder
        self._mask = np.array([i in self.post_select_det_ids 
                                for i in range(self.num_detectors)],dtype=np.uint8)
        if self.decoder is None:
            # self.decoder = pymatching.Matching(self.dem)
            raise ValueError('decoder unset')
        

    @staticmethod
    def from_task(task: sinter.Task) -> CompiledSOXSampler:
        stim_circ = task.circuit
        dem_helper = DEM(stim_circ)
        dem = dem_helper.prune_post_selected()
        post_select_det_ids = frozenset(dem_helper.post_selection_mask)
        decoder = construct_decoder(task, basis='X')
        return CompiledSOXSampler(dem,stim_circ,post_select_det_ids,decoder)
    

    def sample(self, suggested_shots) -> sinter.AnonTaskStats:
        t0 = time.monotonic()
        dets, obs = self.sampler.sample(shots=suggested_shots, separate_observables=True)
        keep_mask = ~np.any(dets & self._mask, axis=1)
        dets = dets[keep_mask]
        obs = np.reshape(obs[keep_mask],(-1))
        predicted_obs, soft_output = self.decoder.decode_batch_soft_output(dets,return_weights=True)
        soft_output_rounded = np.round(np.reshape(soft_output,(-1))*10/(2.303)).astype(np.int64)
        errs = np.reshape(predicted_obs,(-1)) ^ obs
        counter = collections.Counter()
        for gap, err in zip(soft_output_rounded, errs):
            counter[f'E{gap}' if err else f'C{gap}'] += 1
        t1 = time.monotonic()

        return sinter.AnonTaskStats(
            shots=suggested_shots,
            errors=np.count_nonzero(errs),
            discards=suggested_shots - np.count_nonzero(keep_mask),
            seconds=t1 - t0,
            custom_counts=counter,
        )




def construct_decoder(task: sinter.Task, basis: str) -> pymatching.Matching:
    if basis not in ['X','Z','BOTH']:
        raise ValueError('basis should be one of '
                         '{X,Z,BOTH}')
    dem = DEM(task.circuit)

    pruned_dem = dem.prune_post_selected()

    unpost_det_ids = []
    unpost_time = []
    for i in range(dem.num_dets):
        if i not in dem.post_selection_mask:
            unpost_det_ids.append(i)
            unpost_time.append(round(dem.det_id_coords[i][-1]))


    x_upper_bd_dets = []
    x_lower_bd_dets = []
    z_left_bd_dets = []
    z_right_bd_dets = []



    coords_to_det_id = {}
    for det_id in unpost_det_ids:
        det_id_coords = dem.det_id_coords[det_id]
        coords_to_det_id[tuple(det_id_coords)] = det_id
        try:
            err_dets_list = dem.det_id_connection[det_id]
        except:
            continue
        for err_dets in err_dets_list:
            if len(err_dets) == 1:
                if (det_id_coords[0] - det_id_coords[1]) % 2 == 0:
                    if det_id_coords[1] < 0:
                        x_upper_bd_dets.append(det_id)
                        # print(det_id_coords,'xupper')
                        break
                    elif det_id_coords[1] > 0:
                        x_lower_bd_dets.append(det_id)
                        # print(det_id_coords,'xlower')
                        break
                if (det_id_coords[0] - det_id_coords[1]) % 2 == 1:
                    if det_id_coords[0] < 0:
                        z_left_bd_dets.append(det_id)
                        # print(det_id_coords,'zleft')
                        break
                    elif det_id_coords[0] > 0:
                        z_right_bd_dets.append(det_id)
                        # print(det_id_coords,'zright')
                        break

    if basis == 'X':
        x_upper_index = dem.num_dets
        x_lower_index = x_upper_index + 1
    elif basis == 'Z':
        z_left_index == dem.num_dets
        z_right_index = z_left_index + 1
    else:
        assert basis == 'BOTH'
        x_upper_index = dem.num_dets
        x_lower_index = x_upper_index + 1
        z_left_index = x_lower_index + 1
        z_right_index = z_left_index + 1



    decoder = pymatching.Matching(pruned_dem)
    decoder.SO_calculator_setup()
    if basis == 'X':
        for bd_index in [x_upper_index,x_lower_index]:
            decoder.add_boundary_node_SO(bd_index)
    elif basis == 'Z':
        for bd_index in [z_left_index,z_right_index]:
            decoder.add_boundary_node_SO(bd_index)
    else:
        assert basis == 'BOTH'
        for bd_index in [x_upper_index,x_lower_index,
                         z_left_index,z_right_index]:
            decoder.add_boundary_node_SO(bd_index)
    

    if basis == 'Z' or basis == 'BOTH':
        for bd_det_id in z_left_bd_dets:
            decoder.add_boundary_edge_SO(bd_det_id,z_left_index)
        for bd_det_id in z_right_bd_dets:
            decoder.add_boundary_edge_SO(bd_det_id,z_right_index)
    if basis == 'X' or basis == 'BOTH':
        for bd_det_id in x_lower_bd_dets:
            decoder.add_boundary_edge_SO(bd_det_id,x_lower_index)
        for bd_det_id in x_upper_bd_dets:
            decoder.add_boundary_edge_SO(bd_det_id,x_upper_index)

    if basis == 'Z' or basis == 'BOTH':
        decoder.add_cycle_endpoints_pair_SO(z_left_index,z_right_index)
    if basis == 'X' or basis == 'BOTH':
        decoder.add_cycle_endpoints_pair_SO(x_lower_index,x_upper_index)

    return decoder            