import vpython_noGUI as vpy
#import vpython as vpy
import numpy as np
import interaction_energy
import Molecules

def create_random_origin(side, radius_of_disk, no_of_discretization, buffer):
    
    good_zone_lower_limit = []
    good_zone_upper_limit = []
    
    step = side/no_of_discretization
    
    x_pos = 0
    
    for i in range(no_of_discretization):
    
        #print(x_pos-buffer, x_pos+buffer)
        good_zone_lower_limit.append(x_pos-buffer)
        good_zone_upper_limit.append(x_pos+buffer)
        x_pos += step
        if x_pos > side: x_pos -= side
    
    #print("Limits")
    
    random_x_array = []
    random_y_array = []
    
    for i in range(int(len(good_zone_upper_limit)-1)):
        
        #print("Acessing ", i, i+1, "out of ", len(good_zone_upper_limit))
        
        # else:
        #print(good_zone_upper_limit[i], " to ", good_zone_lower_limit[i+1])
            
        random_x, random_y  = np.random.default_rng().uniform(low=good_zone_upper_limit[i],
                                                              high=good_zone_lower_limit[i+1], size=2)
        
        random_x_array.append(random_x)
        random_y_array.append(random_y)
        
        if i == len(good_zone_upper_limit)-2:         
            
            #print(good_zone_upper_limit[i+1], " to ", side-buffer)

            random_x, random_y  = np.random.default_rng().uniform(low=good_zone_upper_limit[i+1],
                                                                  high=side-buffer, size=2)
            
            random_x_array.append(random_x)
            random_y_array.append(random_y)
        
    #print("")
    random_x_index, random_y_index = np.random.choice(len(random_x_array), size=2)
    
    x_pos = random_x_array[random_x_index]
    y_pos = random_y_array[random_y_index]
    
    random_origin = vpy.sphere(pos=vpy.vec(x_pos, y_pos, 0), radius=4*radius_of_disk, color=vpy.color.magenta, visible=False)
    
    # for i in range(no_of_discretization):
    #     for j in range(no_of_discretization):
            
    #         print(x_pos, y_pos)
        
    #         grid_line_x = vpy.curve(vpy.vec(x_pos,0,0), vpy.vec(x_pos,side,0), radius=0, color=vpy.color.black)
    #         grid_line_y = vpy.curve(vpy.vec(0,y_pos,0), vpy.vec(side,y_pos,0), radius=0, color=vpy.color.black)
            
    
    #         vpy.sphere(pos=vpy.vec(x_pos,y_pos,0), radius=radius_of_disk, color=vpy.color.black)
            
    #         #x_upper = x_pos+buffer
            
    #         vpy.sphere(pos=vpy.vec(x_pos+buffer,y_pos+buffer,0), radius=radius_of_disk, color=vpy.color.green)
    #         vpy.sphere(pos=vpy.vec(x_pos-buffer,y_pos+buffer,0), radius=radius_of_disk, color=vpy.color.green)
    #         #x_lower = x_pos-buffer
        
            
    #         vpy.sphere(pos=vpy.vec(x_pos+buffer,y_pos-buffer,0), radius=radius_of_disk, color=vpy.color.green)
    #         vpy.sphere(pos=vpy.vec(x_pos-buffer,y_pos-buffer,0), radius=radius_of_disk, color=vpy.color.green)
            
    #         #print(x_pos)
            
    #         y_pos += step
    #         if y_pos > side: y_pos -= side
            
    #     x_pos += step
    #     if x_pos > side: x_pos -= side
    
    return random_origin
    
def create_grid_from_random_origin(random_origin, side, no_of_discretization, radius_of_disk):
    
    step = side/no_of_discretization
    
    x_pos = random_origin.pos.x
    y_pos = random_origin.pos.y
    
    
    
    cell_dim_dict = {}
    
    for i in range(no_of_discretization):
        
        
        for j in range(no_of_discretization):
            
            #cell_border_dims_array = []
            
            #print(x_pos, y_pos)
        
            #grid_line_x = vpy.curve(vpy.vec(x_pos,0,0), vpy.vec(x_pos,side,0), radius=0, color=vpy.color.black)
            #grid_line_y = vpy.curve(vpy.vec(0,y_pos,0), vpy.vec(side,y_pos,0), radius=0, color=vpy.color.black)
            
            x_lower = x_pos
            y_lower = y_pos

            x_upper = x_pos+step
            y_upper = y_pos+step
            
            x_split_right = (x_lower+x_upper)/2
            y_split_right = (y_lower+y_upper)/2
            
            x_split_left = (x_lower+x_upper)/2
            y_split_left = (y_lower+y_upper)/2
            
            if x_upper > side:
                x_upper -= side
                x_split_right = side
                x_split_left  = 0
                
            if y_upper > side:
                y_upper -= side
                y_split_right = side
                y_split_left  = 0
            
            # vpy.sphere(pos=vpy.vec(x_lower, y_lower, 0), radius=0.5, color=vpy.color.blue)
            # vpy.sphere(pos=vpy.vec(x_lower, y_upper, 0), radius=0.5, color=vpy.color.blue)
            # vpy.sphere(pos=vpy.vec(x_upper, y_lower, 0), radius=0.5, color=vpy.color.red)
            # vpy.sphere(pos=vpy.vec(x_upper, y_upper, 0), radius=0.5, color=vpy.color.red) 
            
            # vpy.sphere(pos=vpy.vec(x_split_left, y_lower, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_split_left, y_upper, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_lower, y_split_left, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_upper, y_split_left, 0), radius=0.5, color=vpy.color.black)
            
            # vpy.sphere(pos=vpy.vec(x_split_right, y_lower, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_split_right, y_upper, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_lower, y_split_right, 0), radius=0.5, color=vpy.color.black)
            # vpy.sphere(pos=vpy.vec(x_upper, y_split_right, 0), radius=0.5, color=vpy.color.black)
            
            #2 zones are defined for each axis x_lower to x_split_right and x_split_left to x_higher
            
            
            #####################################################################################
            
            x_lower_border = x_lower + 2*radius_of_disk
            x_upper_border = x_upper - 2*radius_of_disk
            
            y_lower_border = y_lower + 2*radius_of_disk
            y_upper_border = y_upper - 2*radius_of_disk
                        
            # vpy.sphere(pos=vpy.vec(x_lower_border, y_lower_border, 0), radius=0.5, color=vpy.color.green)
            # vpy.sphere(pos=vpy.vec(x_lower_border, y_upper_border, 0), radius=0.5, color=vpy.color.green)
            # vpy.sphere(pos=vpy.vec(x_upper_border, y_lower_border, 0), radius=0.5, color=vpy.color.green)
            # vpy.sphere(pos=vpy.vec(x_upper_border, y_upper_border, 0), radius=0.5, color=vpy.color.green) 

            
            #print(system_molecule_array, "\n \n", distance_array)
            


            
            assert x_split_right == side or x_split_right == x_split_left
            assert x_split_left == 0 or x_split_left == x_split_right
            
            assert 0 < y_lower_border < side
            assert 0 < y_upper_border < side
            assert 0 < x_lower_border < side
            assert 0 < x_upper_border < side
    
            cell_border_dims_array = [[x_lower, x_upper,
                                       y_lower, y_upper],
                                     [x_split_left, x_split_right,
                                       y_split_left, y_split_right],
                                     [x_lower_border, x_upper_border,
                                       y_lower_border, y_upper_border]]
            
            #print(cell_border_dims_array)            
            #print(x_pos)
            
            y_pos += step
            if y_pos > side: y_pos -= side
        
            cell_dim_dict[str(i)+"_"+str(j)] = cell_border_dims_array
            
        x_pos += step
        if x_pos > side: x_pos -= side
    

    return cell_dim_dict
        

def grid_up_the_system(system_array_xyz, system_array_gui, no_of_discretization, side, radius_of_disk):

    
    step = side/no_of_discretization #8
    buffer = 4*radius_of_disk
    
    
    random_origin = create_random_origin(side, radius_of_disk, no_of_discretization, buffer)
    
    cell_dim_dict = create_grid_from_random_origin(random_origin, side, no_of_discretization, radius_of_disk)
    
    #for key in cell_dim_dict:
    #    print("cell_border_dict", key, cell_dim_dict[key], "\n")
    
    diagonal = np.linalg.norm(np.array([step/2, step/2, 0])-np.array([0, 0, 0]))
    
    cell_center_array = []
    
    cell_dict = {}
    cell_border_dict = {}
    
    # going over every cell
    for i in range(no_of_discretization):
        for j in range(no_of_discretization):
            
            
            cell_border_dims_array = cell_dim_dict[str(i) +"_" +str(j)]
            
            #print("cell_border_dims_array", cell_border_dims_array)
            
            [x_lower, x_upper, y_lower, y_upper] = cell_border_dims_array[0]
            [x_split_left, x_split_right, y_split_left, y_split_right] = cell_border_dims_array[1]
            [x_lower_border, x_upper_border, y_lower_border, y_upper_border] = cell_border_dims_array[2]
    
            #print(i,j)
            
            cell_sphere_array = []
            cell_border_sphere_array = []
                        
            if (i%2 == 0 and j%2 == 0) or (i%2 == 1 and j%2 == 1): color_array = [0,0,0]
            else: color_array = [255,255,255]
            
            color = vpy.vec(color_array[0]/255,color_array[1]/255,color_array[2]/255)
            
            #vpy.sphere(pos=vpy.vec(x_b, y_l, 0), radius=2, color=color)
            
            cell_center_x, cell_center_y = x_lower+step/2, y_lower+step/2
            
            if cell_center_x > side: cell_center_x -= side
            if cell_center_y > side: cell_center_y -= side
            
            cell_center_array.append([cell_center_x, cell_center_y, 0])
            
            #vpy.sphere(pos=vpy.vec(cell_center_x, cell_center_y, 0), radius=1, color=color)
    
            #vpy.sphere(pos=vpy.vec(cell_center_x, cell_center_y, 0), radius=diagonal, opacity=0.05, color=color)
            
            #print(x_b, y_l, "  to  ",  x_b+x_b/2, y_l+y_l/2)
        
            cell_center = np.array([np.array([cell_center_x, cell_center_y, 0])])
            
            neighbour_list = interaction_energy.build_neighbour_list(system_molecules_array=system_array_xyz,
                                                                      trial_molecules_array=cell_center,
                                                                      box_dims=[side, side, side],
                                                                      cut_off=diagonal)
            
            system_molecule_array, trial_molecule_array, distance_array = neighbour_list[0]


            for x,y,z in system_molecule_array:
                

                is_inside_x_limits_bool = x_lower <= x <= x_split_right or x_split_left <= x <= x_upper
    
                is_inside_y_limits_bool = y_lower <= y <= y_split_right or y_split_left <= y <= y_upper
                
                if is_inside_x_limits_bool and is_inside_y_limits_bool:
                    

                    left_bool = x_lower <= x <= x_lower_border
                    right_bool = x_upper_border <= x <= x_upper
                    bottom_bool = y_lower <= y <= y_lower_border
                    top_bool = y_upper_border <= y <= y_upper
                    
                    
                    if  left_bool or right_bool or top_bool or bottom_bool:
                        #simple_sphere = vpy.simple_sphere(pos=vpy.vec(x,y,z), radius=radius_of_disk, color=vpy.color.red)
                        atom = Molecules.Atom(atomic_symbol="D", position=[x,y,z])
                        atom.radius= radius_of_disk
                        atom.color = vpy.color.red
                        cell_border_sphere_array.append(atom)
                    else:    
                        atom = Molecules.Atom(atomic_symbol="D", position=[x,y,z])
                        atom.radius= radius_of_disk
                        atom.color = color
                        cell_sphere_array.append(atom)
                    
                        
            cell_dict[str(i)+"_"+str(j)] = cell_sphere_array
            cell_border_dict[str(i)+"_"+str(j)] = cell_border_sphere_array
            #cell_border_dims[str(i)+"_"+str(j)] = cell_border_dims_array

    #print("Total_spheres = ", len(system_array_xyz))
    #print("border_sphers = ", len(cell_border_dict))
    
    Total_cell_spheres = 0
    for keys in cell_dict: Total_cell_spheres += len(cell_dict[keys])
    #print("Total cell spheres", Total_cell_spheres)


    Total_cell_border_spheres = 0    
    for keys in cell_border_dict: Total_cell_border_spheres += len(cell_border_dict[keys])
    #print("Total cell border spheres", Total_cell_border_spheres)

    #assert len(system_array_xyz) == Total_cell_border_spheres + Total_cell_spheres
    
    return cell_border_dict, cell_dict, cell_dim_dict


        


if __name__ == '__main__':

    screen_width : int = 2500
    screen_height : int = 1000

    vpy.scene.background = vpy.color.white
    vpy.scene.width = screen_width
    vpy.scene.height = screen_height

    #Disables right click spin by user's mouse
    vpy.scene.userspin = False

    no_of_disk = 1000 #107620
    radius_of_disk = 0.5
    side = 50 #350
    
    buffer = 4*radius_of_disk
    
    no_of_discretization = 2
    
    step = side/no_of_discretization
    #Making sure the camera focuses on center of the sheet
    vpy.scene.camera.pos = vpy.vec(side/2,side/2,side/2)
    #Zooming in a bit
    vpy.scene.range = 25

    def create_sheet(side):
        
        corner_bl = vpy.vec(0,0,0)
        corner_br = vpy.vec(0,side,0)
        corner_tl = vpy.vec(side,0,0)
        corner_tr = vpy.vec(side,side,0)
        
        side_bottom = vpy.curve(corner_bl, corner_br, radius=0, color=vpy.color.black)
        side_top    = vpy.curve(corner_tl, corner_tr, radius=0, color=vpy.color.black)
        side_left   = vpy.curve(corner_bl, corner_tl, radius=0, color=vpy.color.black)
        side_right  = vpy.curve(corner_br, corner_tr, radius=0, color=vpy.color.black)
        
        return [side_bottom, side_top, side_left, side_right]
    

    
    create_sheet(side)

    ##########################3

    system_array_xyz = np.random.default_rng().uniform(low=0+radius_of_disk,
                                                          high=side -(radius_of_disk),
                                                          size=[2, no_of_disk])
    system_array_xyz = np.insert(system_array_xyz, 2, 0, axis=0)

    system_array_xyz = system_array_xyz.T

    system_array_gui = []
    for x,y,z in system_array_xyz:
        simple_sphere = vpy.simple_sphere(pos=vpy.vec(x,y,z), radius=radius_of_disk, color=vpy.color.green, opacity=0.1)
        system_array_gui.append(simple_sphere)
        
    #############################
    
    #random_origin = create_random_origin(side, radius_of_disk, step, buffer)
    
    #cell_border_dims = create_grid_from_random_origin(random_origin, step)

    # vpy.curve(vpy.vec(side/2,0,0), vpy.vec(side/2,side,0), radius=0.1, color=vpy.color.black)
    # vpy.curve(vpy.vec(0,side/2,0), vpy.vec(side,side/2,0), radius=0.1, color=vpy.color.black)
    
    # x_lower = 25
    # y_lower = 25
    # x_upper = 75
    # y_upper = 75
        
    # for x,y,z in system_array_xyz:

    #     x,y = apply_peroidic_boundary_conditions(x, y,
    #                                               x_lower, y_lower,
    #                                               x_upper, y_upper,
    #                                               side_x=side, side_y=side)

    #     # x,y = apply_peroidic_boundary_conditions(x, y,
    #     #                                           x_lower=75, y_lower=75,
    #     #                                           x_upper=25, y_upper=25,
    #     #                                           side_x=100, side_y=100)
        
    #     # x,y = apply_peroidic_boundary_conditions(x, y,
    #     #                                           x_lower=12.5, y_lower=12.5,
    #     #                                           x_upper=87.5, y_upper=87.5,
    #     #                                           side_x=100, side_y=100)
        
        
    #     simple_sphere = vpy.simple_sphere(pos=vpy.vec(x,y,z), radius=radius_of_disk, color=vpy.color.red)
    # #draw_cube_from_points(x_start=87.5, x_end=12.5, y_start, y_end, z_start, z_end)
    
    # vpy.sphere(pos=vpy.vec(x_lower,y_lower,0), radius=1, color=vpy.color.blue)
    # vpy.sphere(pos=vpy.vec(x_lower,y_upper,0), radius=1, color=vpy.color.blue)
    # vpy.sphere(pos=vpy.vec(x_upper,y_lower,0), radius=1, color=vpy.color.blue)
    # vpy.sphere(pos=vpy.vec(x_upper,y_upper,0), radius=1, color=vpy.color.blue)    
    
    ################################3
    
    #center = vpy.sphere(pos=vpy.vec(side/2+side/12,side/2+side/12,0), radius=2, color=vpy.color.blue)


    #for i in range(100):    
    cell_border_dict, cell_dict, cell_dim_dict = grid_up_the_system(system_array_xyz, system_array_gui, no_of_discretization, side, radius_of_disk)

    #for key in cell_border_dims:
    #    print(key, cell_border_dims[key])

    #for cell_sphere_array in cell_dict.values(): print(len(ce))         
                
                #simple_sphere = vpy.simple_sphere(pos=vpy.vec(x,y,z), radius=radius_of_disk, color=color)
                        
            #system_array_gui.append(simple_sphere)
        
        
    


