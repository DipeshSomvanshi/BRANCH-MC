import numpy as np
from scipy.spatial import cKDTree
import os
import glob
import argparse
import time

parser = argparse.ArgumentParser(description='GCMC for hard disks.')

parser.add_argument('--restart', action='store_true', help='Restarts the simultion from last saved file')

parser.add_argument('--chem_pot', type=float, help='Chemical potential of the system')

args = parser.parse_args()

restart_mode = args.restart
print(f"restart mode: {restart_mode}")

if args.chem_pot is None:
    parser.error("Please provide a chemical potential with --chem_pot flag.")
else:
    chem_pot = float(args.chem_pot)

# chem_pot = 1
# restart_mode = False

print(f"chem_pot: {chem_pot}")

fugacity = np.exp(chem_pot)
print("fugacity: ", fugacity, "chem_pot: ", chem_pot, "add chances", 1/fugacity)

file_name = 'rp_' + str(chem_pot) + '.out'
print(file_name)

burn_in = 0

if chem_pot > 9:
    mc_cycles = 10.0025e9
else:
    mc_cycles = 1.0025e9

side_of_sheet = 10

no_neigh = 10

move_factor = 0.25

radius_of_disk = 0.5

volume_of_sphere = 4/3 * np.pi * radius_of_disk * radius_of_disk * radius_of_disk
Total_volume = side_of_sheet * side_of_sheet * side_of_sheet

molecules_xyz = np.array([[side_of_sheet/2 , side_of_sheet/2 , side_of_sheet/2]])

box_dims = [side_of_sheet , side_of_sheet , side_of_sheet]

save_file_folder = "run-files"

def run_mc_in_a_cell(molecules_xyz:list, mc_cycles):

    if restart_mode == False:

        print("Starting a new sim")

        add_accept_counter = 0
        add_step_counter = 0

        move_step_counter = 0
        move_accept_counter = 0

        remove_accept_counter = 0
        remove_step_counter = 0

        repack_accept_counter = 0
        repack_step_counter = 0

        time_per_step = 0

        total_step_counter = 0

        density = (len(molecules_xyz)*volume_of_sphere)/Total_volume    
        old_file_path = 'xyz_chem_pot'+str(chem_pot)+"_cycles_"+str( (total_step_counter - burn_in)/1e6) + '_density_' + str(np.round(density,8)) +  "_add_accept_" + str(add_accept_counter) + "_add_step_" + str(add_step_counter) + "_remove_accept_" + str(remove_accept_counter) + "_remove_step_" + str(remove_step_counter) + "_move_accept_" + str(move_accept_counter) + "_move_step_" + str(move_step_counter) + "_repack_accept_" + str(repack_accept_counter) + "_repack_step_" + str(repack_step_counter) + "_time_per_step_" + str(np.round(time_per_step,3)) + '.xyz'

    else:

        all_xyz_files = glob.glob(save_file_folder +"/*.xyz")

        restart_file_name = None
        for file in all_xyz_files:
            print()
            if float(file.split("_")[3]) == chem_pot:
                restart_file_name = file
                break
        
        if restart_file_name is None:
            print("No file found to restart from. Check the chemical potential and try again. Exiting...")
            return
        
        print("Restarting from file: ", restart_file_name, "\n")

        restart_file_split = restart_file_name.split("_")

        molecules_xyz = np.loadtxt(restart_file_name)
        
        total_step_counter  = float(restart_file_split[5])*1e6
        add_accept_counter    = int(restart_file_split[10])
        add_step_counter      = int(restart_file_split[13])
        move_step_counter     = int(restart_file_split[22])
        move_accept_counter   = int(restart_file_split[25])
        remove_accept_counter = int(restart_file_split[16])
        remove_step_counter   = int(restart_file_split[19])
        repack_accept_counter = int(restart_file_split[28])
        repack_step_counter   = int(restart_file_split[31])
        time_per_step       = float(restart_file_split[35])
        
        density = (len(molecules_xyz)*volume_of_sphere)/Total_volume
        old_file_path = restart_file_name

        print("MC cycle: ", total_step_counter/1e6, "density: ", np.round(density,4), "add_accept: ", int(add_accept_counter), "add_step: ", int(add_step_counter), "remove_accept: ", int(remove_accept_counter), "remove_step: ", int(remove_step_counter), "repack_accept: ", int(repack_accept_counter), "repack_step: ", int(repack_step_counter), "move_accept: ", int(move_accept_counter), "move_step: ", int(move_step_counter), 'time_per_step:', np.round(time_per_step,3), "\n")
    
    start_time = time.time()

    while total_step_counter < mc_cycles:

            
        choice = np.random.choice(["move", "add", "remove"],
                                1,
                                p=[0.5, 0.25, 0.25],
                                replace=True)[0]

        ######### Cavity ###########            

        #Draw a random point in the box
        # random_point = np.random.default_rng().uniform(low=0,
        #                                         high=side_of_sheet,
        #                                         size=3)

        ################## Cavity
        if choice == "move" and len(molecules_xyz) > 0:

            if total_step_counter >= burn_in: move_step_counter += 1

            #Randomly choose a molecule to move
            random_disk_index = np.random.choice(len(molecules_xyz))
            random_disk_xyz = np.array(molecules_xyz[random_disk_index])
            
            #Splitting the system molecules in two parts:
            #Random disk and all sys mol except random disk
            disk_array_without_random_disk_xyz = np.delete(molecules_xyz, random_disk_index, axis=0)

            assert len(disk_array_without_random_disk_xyz) == len(molecules_xyz)-1, str(len(disk_array_without_random_disk_xyz)) + " " + str(len(molecules_xyz))
            
            [ranf_1, ranf_2, ranf_3] = np.random.uniform(low=0, high=1, size=(3))

            moved_disk_x = random_disk_xyz[0]+move_factor*(ranf_1-0.5)
            moved_disk_y = random_disk_xyz[1]+move_factor*(ranf_2-0.5)
            moved_disk_z = random_disk_xyz[2]+move_factor*(ranf_3-0.5)

            if moved_disk_x < 0: moved_disk_x += side_of_sheet
            if moved_disk_x > side_of_sheet: moved_disk_x -= side_of_sheet
            
            if moved_disk_y < 0: moved_disk_y += side_of_sheet
            if moved_disk_y > side_of_sheet: moved_disk_y -= side_of_sheet

            if moved_disk_z < 0: moved_disk_z += side_of_sheet
            if moved_disk_z > side_of_sheet: moved_disk_z -= side_of_sheet

            moved_disk_xyz = np.array([moved_disk_x, moved_disk_y, moved_disk_z])            

            assert 0 < moved_disk_xyz[0] < box_dims[0], "0 < " + str(moved_disk_xyz[0]) + " < " + str(box_dims[0])
            assert 0 < moved_disk_xyz[1] < box_dims[1], "0 < " + str(moved_disk_xyz[1]) + " < " + str(box_dims[1])
            assert 0 < moved_disk_xyz[2] < box_dims[2], "0 < " + str(moved_disk_xyz[2]) + " < " + str(box_dims[2])

            # Construct the KD-tree with periodic boundary conditions
            kdtree = cKDTree(disk_array_without_random_disk_xyz, boxsize=box_dims)

            distances_raw, _ = kdtree.query([moved_disk_xyz],
                                            k=no_neigh,
                                            distance_upper_bound=1.1)

            distances = distances_raw[0]

            distances = distances[distances <= 1.0]

            if len(distances) == 0:

                # if len(molecules_xyz) > 1:

                #     kdtree_cav = cKDTree(disk_array_without_random_disk_xyz, boxsize=box_dims)

                #     distances_raw_cav, _ = kdtree_cav.query([moved_disk_xyz],
                #                                     k=no_neigh,
                #                                     distance_upper_bound=1.1)

                #     distances_cav = distances_raw_cav[0]

                #     distances_cav = distances_cav[distances_cav <= 1.0]
                    
                #     if len(distances_cav) != 0: continue
                
                old_n = len(molecules_xyz)

                molecules_xyz = np.append(disk_array_without_random_disk_xyz, [moved_disk_xyz], axis=0)
                
                assert len(molecules_xyz) == old_n, str(len(molecules_xyz)) + " " + str(old_n)

                if total_step_counter >= burn_in: move_accept_counter += 1

        elif choice == "add":

            if total_step_counter >= burn_in: add_step_counter += 1
            total_step_counter += 1

            if total_step_counter % 5e3 == 0:
                
                time_per_step = time.time() - start_time

                density = (len(molecules_xyz)*volume_of_sphere)/Total_volume            
                print("MC cycle: ", total_step_counter/1e6, "density: ", np.round(density,4), "add_accept: ", int(add_accept_counter), "add_step: ", int(add_step_counter), "remove_accept: ", int(remove_accept_counter), "remove_step: ", int(remove_step_counter), "repack_accept: ", int(repack_accept_counter), "repack_step: ", int(repack_step_counter), "move_accept: ", int(move_accept_counter), "move_step: ", int(move_step_counter), 'time_per_step:', np.round(time_per_step,3), "\n")

                if total_step_counter >= burn_in:
                    with open(file_name, 'a') as f:
                        f.write(str( (total_step_counter - burn_in)/1e6) + '_density_' + str(np.round(density,8)) +  "_add_accept_" + str(add_accept_counter) + "_add_step_" + str(add_step_counter) + "_remove_accept_" + str(remove_accept_counter) + "_remove_step_" + str(remove_step_counter) + "_move_accept_" + str(move_accept_counter) + "_move_step_" + str(move_step_counter) + "_repack_accept_" + str(repack_accept_counter) + "_repack_step_" + str(repack_step_counter) + "_time_per_step_" + str(np.round(time_per_step,3)) + '\n')
                
                if os.path.exists(old_file_path):os.remove(old_file_path)

                density = (len(molecules_xyz)*volume_of_sphere)/Total_volume
                
                new_file_name = 'xyz_chem_pot_'+str(chem_pot)+"_cycles_"+str( (total_step_counter - burn_in)/1e6) + '_density_' + str(np.round(density,8)) +  "_add_accept_" + str(add_accept_counter) + "_add_step_" + str(add_step_counter) + "_remove_accept_" + str(remove_accept_counter) + "_remove_step_" + str(remove_step_counter) + "_move_accept_" + str(move_accept_counter) + "_move_step_" + str(move_step_counter) + "_repack_accept_" + str(repack_accept_counter) + "_repack_step_" + str(repack_step_counter) + "_time_per_step_" + str(np.round(time_per_step,3)) + '_.xyz'

                # Save the new xyz file
                new_file_path = os.path.join(save_file_folder, new_file_name)
                np.savetxt(new_file_path, molecules_xyz, fmt='%s')
                
                old_file_path = new_file_path

                start_time = time.time()

            # trial_molecule_x, trial_molecule_y = np.random.uniform(low=0, high=side_of_sheet, size=(1, 2))[0]
            # trial_molecule_z = 0

            ################ Generate random molecule in the box ################

            [trial_molecule_x, trial_molecule_y, trial_molecule_z] = np.random.uniform(0, side_of_sheet, 3)

            if trial_molecule_x <= 0: trial_molecule_x += side_of_sheet
            if trial_molecule_x > side_of_sheet: trial_molecule_x -= side_of_sheet
            
            if trial_molecule_y <= 0: trial_molecule_y += side_of_sheet
            if trial_molecule_y > side_of_sheet: trial_molecule_y -= side_of_sheet

            if trial_molecule_z <= 0: trial_molecule_z += side_of_sheet
            if trial_molecule_z > side_of_sheet: trial_molecule_z -= side_of_sheet

            trial_molecule_xyz = np.array([trial_molecule_x, trial_molecule_y, trial_molecule_z])

            ################ Check for overlap ################

            kdtree = cKDTree(molecules_xyz, boxsize=box_dims)

            distances_raw, _ = kdtree.query([trial_molecule_xyz],
                                                       k=no_neigh,
                                                       distance_upper_bound=1.1)

            distances = distances_raw[0]

            distances = distances[distances <= 1.0]

            if len(distances) == 0:
                
                alpha = Total_volume/(len(molecules_xyz)+1) * fugacity

                random_number = np.random.uniform(low=0, high=1, size=(1))[0]
                if alpha > random_number:

                    molecules_xyz = np.append(molecules_xyz, [trial_molecule_xyz], axis=0)
                    if total_step_counter >= burn_in: add_accept_counter += 1
        
        elif choice == "remove":

            if len(molecules_xyz) == 1: continue

            if total_step_counter >= burn_in: remove_step_counter += 1
            
            total_step_counter += 1

            if total_step_counter % 5e3 == 0:
                
                time_per_step = time.time() - start_time

                density = (len(molecules_xyz)*volume_of_sphere)/Total_volume            
                print("MC cycle: ", total_step_counter/1e6, "density: ", np.round(density,4), "add_accept: ", int(add_accept_counter), "add_step: ", int(add_step_counter), "remove_accept: ", int(remove_accept_counter), "remove_step: ", int(remove_step_counter), "repack_accept: ", int(repack_accept_counter), "repack_step: ", int(repack_step_counter), "move_accept: ", int(move_accept_counter), "move_step: ", int(move_step_counter), 'time_per_step:', np.round(time_per_step,3), "\n")

                if total_step_counter >= burn_in:
                    with open(file_name, 'a') as f:
                        f.write(str( (total_step_counter - burn_in)/1e6) + '_density_' + str(np.round(density,8)) +  "_add_accept_" + str(add_accept_counter) + "_add_step_" + str(add_step_counter) + "_remove_accept_" + str(remove_accept_counter) + "_remove_step_" + str(remove_step_counter) + "_move_accept_" + str(move_accept_counter) + "_move_step_" + str(move_step_counter) + "_repack_accept_" + str(repack_accept_counter) + "_repack_step_" + str(repack_step_counter) + "_time_per_step_" + str(np.round(time_per_step,3)) + '\n')
                
                if os.path.exists(old_file_path):os.remove(old_file_path)
                
                density = (len(molecules_xyz)*volume_of_sphere)/Total_volume
                new_file_name = 'xyz_chem_pot_'+str(chem_pot)+"_cycles_"+str( (total_step_counter - burn_in)/1e6) + '_density_' + str(np.round(density,8)) +  "_add_accept_" + str(add_accept_counter) + "_add_step_" + str(add_step_counter) + "_remove_accept_" + str(remove_accept_counter) + "_remove_step_" + str(remove_step_counter) + "_move_accept_" + str(move_accept_counter) + "_move_step_" + str(move_step_counter) + "_repack_accept_" + str(repack_accept_counter) + "_repack_step_" + str(repack_step_counter) + "_time_per_step_" + str(np.round(time_per_step,3)) + '_.xyz'
                
                # Save the new xyz file
                new_file_path = os.path.join(save_file_folder, new_file_name)

                np.savetxt(new_file_path, molecules_xyz, fmt='%s')
                
                old_file_path = new_file_path

                start_time = time.time()

            alpha = len(molecules_xyz)/Total_volume * 1/fugacity

            random_number = np.random.uniform(low=0, high=1, size=(1))[0]
            
            if alpha > random_number: 

                random_disk_index = np.random.choice(len(molecules_xyz))

                molecules_xyz = np.delete(molecules_xyz, random_disk_index, axis=0)

                if total_step_counter >= burn_in: remove_accept_counter += 1

run_mc_in_a_cell(molecules_xyz, mc_cycles)
