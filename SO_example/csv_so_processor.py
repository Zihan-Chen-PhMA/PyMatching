import sinter
import numpy as np
import os.path as path


class RES_Handler():

    def __init__(self, 
                 file_name: str):
        self.filename_save = file_name
        self.shots = 0
        self.discards = 0
        self.errors = 0
        pass


    def rate(self, hits_list: list[int], shots_list: list[int]) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=1000).best)
        return ret_list

    def rate_ceil_binfit(self, hits_list: list[int],
                         shots_list: list[int],
                         param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).high)
        return ret_list

    def rate_floor_binfit(self, hits_list: list[int],
                          shots_list: list[int],
                          param: int = 1000) -> list[float]:
        if len(hits_list) != len(shots_list):
            raise ValueError('input lengths mismatch')
        ret_list = []
        for hits, shots in zip(hits_list, shots_list):
            ret_list.append(sinter.fit_binomial(num_shots=shots,
                                num_hits=hits,
                                max_likelihood_factor=param).low)
        return ret_list
    
    def read_through_custom_counts(self) -> None:
        if path.isfile(self.filename_save) == False:
            raise ValueError('no result found')
        self.stats = sinter.read_stats_from_csv_files(self.filename_save)
        if len(self.stats) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        self.stats = self.stats[0]
        self.shots = self.stats.shots
        self.discards = self.stats.discards
        self.errors = self.stats.errors
    

        


class SoftOutput(RES_Handler):

    def __init__(self, file_name: str):
        super().__init__(file_name)
        self.gap_vals = []
        self.at_gap_shots = []
        self.at_gap_hits = []
        self.geq_gap_hits = []
        self.geq_gap_shots = []
        pass

    def read_through_custom_counts(self):
        if path.isfile(self.filename_save) == False:
            raise ValueError('no result found')
        self.stats = sinter.read_stats_from_csv_files(self.filename_save)
        if len(self.stats) != 1:
            raise ValueError('More than one task result is stored. Unsupported behavior for now.')
        self.stats = self.stats[0]
        self.shots = self.stats.shots
        self.discards = self.stats.discards
        self.errors = self.stats.errors
        custom = self.stats.custom_counts
        gap_val_set = set()
        for key in custom:
            gap_val_set.add(int(key[1:]))
        max_gap = max(gap_val_set)
        self.gap_vals = [i for i in range(max_gap+1)]
        for gap in self.gap_vals:
            try:
                custom['C'+str(gap)]
            except:
                custom['C'+str(gap)] = 0
            try:
                custom['E'+str(gap)]
            except:
                custom['E'+str(gap)] = 0
        e_counts = 0
        c_counts = 0
        for gap in reversed(self.gap_vals):
            e_counts += custom['E'+str(gap)]
            c_counts += custom['C'+str(gap)]
            self.at_gap_shots.insert(0,custom['E'+str(gap)]
                                        + custom['C'+str(gap)])
            self.at_gap_hits.insert(0,custom['E'+str(gap)])
            self.geq_gap_hits.insert(0,e_counts)
            self.geq_gap_shots.insert(0,e_counts+c_counts)




