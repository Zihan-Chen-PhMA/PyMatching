from circuit_gadget import *
from so_sampler import *
from csv_so_processor import *


def main():
    d_sc = 5
    sc_name = 'sc' + str(d_sc)
    sc_code = Rotated_surface_code(d_sc, sc_name)
    n_rounds = d_sc*2
    circ_name =  'sc' + str(d_sc) + '_' + str(n_rounds) + '_SE_rds'  + '.stim'
    result_file = 'sc' + str(d_sc) + '_' + str(n_rounds) + '_SE_rds_SO.csv'
    circ_svg_name = 'sc' + str(d_sc) + '_' + str(n_rounds) + '_SE_rds'  + '.svg'
    sc_code.logic_x_selection(0)
    sc_code.logic_z_selection(0)


    noise = 0.001

    cleaness = False
    detector_val = False

    ft_circuit = Meta_Circuit([sc_code],circ_name,noise)

    ft_circuit.X_init(sc_name,cleaness)

    for i in range(n_rounds):
        ft_circuit.SE_round(sc_name, cleaness=cleaness)
    

    ft_circuit.X_meas(sc_name,cleaness)


    stim_circ = ft_circuit.get_stim_circuit()

    with open(circ_svg_name,'w') as file:
        out = stim_circ.diagram('detslice-with-ops-svg')
        file.write(out.__str__())


    max_shots = 1000_000_000

    sinter_task = sinter.Task(circuit=stim_circ)

    samples = sinter.collect(
        num_workers=12,
        max_shots=max_shots,
        # max_errors=1000,
        tasks=[sinter_task],
        decoders=['SOXSampler'],
        custom_decoders={'SOXSampler':SOXSampler()},
        save_resume_filepath=result_file,
        print_progress=True
    )
        



# NOTE: This is actually necessary! If the code inside 'main()' was at the
# module level, the multiprocessing children spawned by sinter.collect would
# also attempt to run that code.
if __name__ == '__main__':
    main()


