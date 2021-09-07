# SPH taichi implementation by mzhang
from pyevtk.hl import pointsToVTK
import taichi as ti
from wave_sph_solver import *
import argparse

# Default run on CPU
# cuda performance has not been tested
ti.init(arch=ti.cpu)


def main():
    dynamic_allocate = False
    save_frames = False
    adaptive_time_step = False
    #method_name = 'WCSPH'
    method_name = 'PCISPH'
    #method_name = 'DFSPH'

    sim_physical_time = 10.0
    max_frame = 30000

    res = (500, 500)
    screen_to_world_ratio = 35
    dx = 0.1
    u, b, l, r = np.array([res[1], 0, 0, res[0]]) / screen_to_world_ratio

    gui = ti.GUI('SPH', res, background_color=0x112F41)
    sph = SPHSolver(res,
                    screen_to_world_ratio, [u, b, l, r],
                    alpha=0.30,
                    dx=dx,
                    max_time=sim_physical_time,
                    max_steps=max_frame,
                    dynamic_allocate=dynamic_allocate,
                    adaptive_time_step=adaptive_time_step,
                    method=SPHSolver.methods[method_name])

    print("Method use: %s" % method_name)
    # Add fluid particles
    #sph.add_cube(lower_corner=[res[0] / 2 / screen_to_world_ratio - 3, 2 * dx],cube_size=[6, 6],velocity=[0.0, 0.0],density=[1000],color=0x068587,material=SPHSolver.material_fluid)
    #sph.add_cube(lower_corner=[0.2, 0.2],cube_size=[11, 11],velocity=[0.0, 0.0],density=[1000],color=0x068587,material=SPHSolver.material_fluid)
    sph.add_file(10000,material=SPHSolver.material_fluid,color=0x068587)
    #print(sph.particle_pressure)
    #print(sph.particle_density)
    #print(sph.particle_velocity)
    #print(sph.particle_positions)


    colors = np.array([0xED553B, 0x068587, 0xEEEEF0, 0xFFFF00],
                      dtype=np.uint32)
    add_cnt = 0.0
    add = False
    save_cnt = 0.0
    output_fps = 60
    save_point = 1.0 / output_fps
    t = 0.0
    frame = 0
    total_start = time.process_time()
    while frame < max_frame and t < sim_physical_time:
        dt = sph.step(frame, t, total_start)
        #print(sph.particle_density.shape)
        particles = sph.particle_info()
        if(frame%10000==0):
            f_velocity="./velocity_%d.txt" % frame
            f_density = "./density_%d.txt" % frame
            f_position="./position_%d.txt" % frame
            f_pressure = "./pressure_%d.txt" % frame
            f_wave_pressure="./wave_pressure_%d.txt" % frame
            arr_velocity=sph.particle_velocity.to_numpy()
            arr_density = sph.particle_density.to_numpy()
            arr_position = sph.particle_positions.to_numpy()
            arr_pressure=sph.particle_pressure.to_numpy()
            arr_wave_preesure=sph.wave_pressure.to_numpy()
            np.savetxt(f_velocity, arr_velocity)
            np.savetxt(f_density, arr_density)
            np.savetxt(f_position, arr_position)
            np.savetxt(f_pressure, arr_pressure)
            np.savetxt(f_wave_pressure,arr_wave_preesure)

        #f_velocity="./velocity/velocity_%d.txt" % frame
        #f_density = "./density/density_%d.txt" % frame
        #f_position="./position/position_%d.txt" % frame
         #f_pressure = "./pressure/pressure_%d.txt" % frame

        #frame1=77777777
        #f_velocity="./velocity/velocity_%d.txt" % frame1
        #f_density = "./density/density_%d.txt" % frame1
        #f_position="./position/position_%d.txt" % frame1
        #f_pressure = "./pressure/pressure_%d.txt" % frame1

        #arr_velocity=sph.particle_velocity.to_numpy()
        #arr_density = sph.particle_density.to_numpy()
        #arr_position = sph.particle_positions.to_numpy()
        #arr_pressure=sph.particle_pressure.to_numpy()
        #np.savetxt(f_velocity, arr_velocity)
        #np.savetxt(f_density, arr_density)
        #np.savetxt(f_position, arr_position)
        #np.savetxt(f_pressure, arr_pressure)

        #print(sph.particle_density)
        #print(1)
        
        



        # if frame == 1000:
        if add:
            sph.add_cube(lower_corner=[4, 9],
                         cube_size=[8.0, 8.0],
                         velocity=[0.0, -10.0],
                         density=[1000.0],
                         color=0xED553B,
                         material=SPHSolver.material_fluid)
        add = False
        if add:
            sph.add_cube(lower_corner=[10, 4],
                         cube_size=[1.0, 1.0],
                         velocity=[0.0, -5.0],
                         density=[1000.0],
                         color=0xED553B,
                         material=SPHSolver.material_fluid)

        # if frame == 1000:
        if add:
            sph.add_cube(lower_corner=[16, 4],
                         cube_size=[2.0, 2.0],
                         velocity=[0.0, -10.0],
                         density=[1000.0],
                         color=0xEEEEF0,
                         material=SPHSolver.material_fluid)
            

        for pos in particles['position']:
            for j in range(len(res)):
                pos[j] *= screen_to_world_ratio / res[j]

        gui.circles(particles['position'],
                    radius=1.5,
                    color=particles['color'])

        if frame % 100 == 0:
            if 1:
                mat = particles["material"]
                fluid = [i for i,m in enumerate(mat) if m]
                
                pos = particles["position"][np.array(fluid)]
                vel = particles["velocity"][np.array(fluid)]
                
                density = particles["density"][np.array(fluid)]
                ddensity = particles["d_density"][np.array(fluid)]
                pressure = particles["pressure"][np.array(fluid)]
                wave_pressure = particles["wave_pressure"][np.array(fluid)]
                pos_x = copyData(pos[:,0])[np.array(fluid)]
                pos_y = copyData(pos[:,1])[np.array(fluid)]
                #pos_z = copyData(pos[:,2])[np.array(fluid)]
                pos_z=np.ones(density.size,dtype=float)
                vel_x = copyData(vel[:,0])[np.array(fluid)]
                vel_y = copyData(vel[:,1])[np.array(fluid)]
                #vel_z = copyData(vel[:,2])[np.array(fluid)]
                vel_z =np.ones(density.size,dtype=float)
                vel = np.linalg.norm(vel, axis=1)[np.array(fluid)]               
    
                pointsToVTK(f'./vtkData/frame_{frame:06d}', pos_x, pos_y, pos_z,
                            data={"vel_x": vel_x, "vel_y": vel_y, "vel_z": vel_z,"vel": vel, "material": mat,
                                  "density": density, "ddensity": ddensity, "pressure": pressure,"wave_pressure":wave_pressure})
        

        # Save in fixed frame interval, for fixed time step
        if not adaptive_time_step and frame % 1 == 0:
            gui.show(f'./pic/{frame:06d}.png' if save_frames else None)

        # Save in fixed frame per second, for adaptive time step
        if adaptive_time_step and save_cnt >= save_point:
            gui.show(f'./pic/{frame:06d}.png' if save_frames else None)
            save_cnt = 0.0

        frame += 1
        t += dt
        save_cnt += dt
        add_cnt += dt

    print('done')

def copyData(value):
    num = value.size
    rvalue = np.zeros(num)
    for i in range(num):
        rvalue[i] = value[i]
    return rvalue

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--method",
                        type=str,
                        default="PCISPH",
                        help="SPH methods: WCSPH, PCISPH, DFSPH")
    opt = parser.parse_args()
    main()
