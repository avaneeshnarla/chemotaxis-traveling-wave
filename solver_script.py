from __future__ import print_function
from mshr import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import ticker
import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
import os
import shutil
import scipy.io
from datetime import datetime

parameters["allow_extrapolation"] = True

upupfolder='big_set_chi_lmax_Feb2021/'
os.makedirs(upupfolder, exist_ok=True)
shutil.copy(__file__, upupfolder+'solver_script.py')
for rho_c in [1e9]:
    upupufolder=upupfolder+'rhoc='+'{0:3.0f}'.format(rho_c)+'/'
    for scaler1 in [1]:#range(10,100,20):
        phi=5
        upfolder=upupufolder+'phi='+str(phi).zfill(3)+'/'
        for scaler in range(100,1001,100):  
            # Units are seconds, OD (for bacteria), mM (for nutrient and attractant), mm
            
            T = 120 *60*60           # final time
            num_steps = (T/2.0)    # number of time steps
            num_steps = int(num_steps)
            dt = T / num_steps # time step size
            lmax=scaler/9/(60*60)
            am=1e-3
            Da=800*1e-6
            D=50*1e-6
            a_0=0.1
            rho0=0.029
            chi_0=D*(phi+1)
            rhoc=rho_c
            Kma=1e-3
            alpha=0.77/(60*60)
            ahill=1
            now = datetime.now()
            current_time = now.strftime("%y_%m_%d_%H_%M_%S")
            folder=upfolder+'lmax'+'{0:4.2e}'.format(lmax)+'_phi_'+str(phi).zfill(4)+'_resln_2000_'+current_time +'/'   
            os.makedirs(folder, exist_ok=True)
            shutil.copy(__file__, folder+'solver_script.py')
            
            # Read mesh from file
            R_max=100.0
            R_min=0.0
            resln=5000
            fn_poly=3;
            mesh = IntervalMesh(resln,R_min,R_max)
            fig=plot(mesh, title='Mesh') 
            dpi=100
            plt.savefig(folder+'mesh.png', dpi=dpi)
            plt.close()
            
            f = open( folder+'parameters.txt', 'w' )
            f.write( 'T = '     + str(T) + '\n' )
            f.write( 'dt = '    + str(dt) + '\n' )
            f.write( 'lmax = '  + str(lmax) + '\n' )
            f.write( 'am = '     + str(am) + '\n' )
            f.write( 'Da = '    + str(Da) + '\n' )
            f.write( 'D = '     + str(D) + '\n' )
            f.write( 'a_0 = '   + str(a_0) + '\n' )
            f.write( 'rho0 = '  + str(rho0) + '\n' )
            f.write( 'chi_0 = '    + str(chi_0) + '\n' )
            f.write( 'rhoc = '     + str(rhoc) + '\n' )
            f.write( 'Kma = '     + str(Kma) + '\n' )
            f.write( 'alpha = '    + str(alpha) + '\n' )
            f.write( 'ahill = '    + str(ahill) + '\n' )
            f.write( 'grid_size = '    + str(R_max-R_min) + '\n' )
            f.write( 'resln = '    + str(resln) + '\n' )
            f.write( 'dpi = '    + str(dpi) + '\n' )
            f.write( 'fn_poly = '    + str(fn_poly) + '\n' )
            f.close()
            
            # Define function space for velocity
            #W = VectorFunctionSpace(mesh, 'P', 3)
            Q = FunctionSpace(mesh, 'P', fn_poly)
            
            # Define function space for system of concentrations
            P1 = FiniteElement('P', interval, fn_poly)
            element = MixedElement([P1, P1])
            V = FunctionSpace(mesh, element)
            
            # Define test functions
            rho_test, a_test = TestFunctions(V)
            
            # Define functions for velocity and concentrations
            w = Function(Q)
            u = Function(V)
            u_n = Function(V)
            
            # Split system functions to access components
            rho_trial, a_trial = split(u)
            rho_iter, a_iter = split(u_n)
            
            u_n.interpolate(Expression(('(tanh((1-x[0]*x[0]*x[0]*x[0]*x[0]*x[0]))+1)*.029/2',
                                        'a0') ,degree=3, a0=a_0, chi=chi_0, D=D, N=10))
            
            # Define source terms
            f_1 = Constant(0.0)
            k = Constant(dt)
            
            # Define expressions used in variational forms
            lam=lmax*(1-rho_trial/rhoc)#*(n_trial/(n_trial+Kmn))#
            mu=alpha*a_trial**ahill/(a_trial**ahill+Kma**ahill)#(alpha1*lam+alpha2)*a_trial/(a_trial+Kma)
            v=(1+a_trial/am)
            w=(chi_0/v*grad(v))
            F = (((rho_trial - rho_iter) / k)*rho_test*dx + div(w*rho_trial)*rho_test*dx 
                + D*dot(grad(rho_trial), grad(rho_test))*dx - lam*rho_trial*rho_test*dx  
                #+ ((n_trial - n_iter) / k)*n_test*dx 
                #+ Dn*dot(grad(n_trial), grad(n_test))*dx + lam*rho_trial*n_test/Y*dx  
                + ((a_trial - a_iter) / k)*a_test*dx
                + Da*dot(grad(a_trial), grad(a_test))*dx 
                + mu*rho_trial*a_test*dx
                - f_1*rho_test*dx - f_1*a_test*dx)
            
            def Left(x, Left):
                return near(x[0],R_min)
            def Right(x, Right):
                return near(x[0],R_max)
    
            #bcrho1 = DirichletBC(V.sub(0), Constant(rhoc), Left)
            bcrho2 = DirichletBC(V.sub(0), Constant(0.0), Right)
            #bca1 = DirichletBC(V.sub(1), Constant(0), Left)
            bca2 = DirichletBC(V.sub(1), Constant(a_0), Right)
            bc=[bcrho2, bca2]
        
            # Time-stepping
            t = 0
                
            peak_location = [0 for y in range(num_steps)]
            max_location = [0 for y in range(num_steps)]
            peak_width = [0 for y in range(num_steps)]
            peak_height = [0 for y in range(num_steps)]
            peak_velocity   = [0 for y in range(num_steps)]
            max_velocity   = [0 for y in range(num_steps)]
            times=[(y*dt)/(60.0) for y in range(num_steps)]
            max_vel_loc=[0 for y in range(num_steps)]
            vstar = [0 for y in range(num_steps)]
            xstar = [0 for y in range(num_steps)]
            rhostar = [0 for y in range(num_steps)]
            dip_location = [0 for y in range(num_steps)]
            dip_velocity = [0 for y in range(num_steps)]
            a_dip = [0 for y in range(num_steps)]
            rho_dip = [0 for y in range(num_steps)]
            for n in range(num_steps):
                times[n]=t;
            
                # Solve variational problem for time step
                try:
                    goal=Function(V)
                    goal.assign(u_n)
                    m=goal.sub(0)*dx()        
                    solve(F == 0, u, bc, tol=1.0e-6, M=m)
                    
                except:
                    break
                
                # Update current time
                t += dt
                print(t)  
                
                _rho_trial, _a_trial = u.split()
                
                vel=project(chi_0/(_a_trial+am)*_a_trial.dx(0))
                rho_g=project(_rho_trial.dx(0))
                arr=np.flip(np.array(rho_g.vector()),0)
                indx=1
                start_indx=0
                peak_indx=1
                for arr_el in range(10,arr.size-10):
                    if (arr[arr_el]<0 and arr[arr_el+1]>0 and start_indx==0):
                        start_indx=(arr_el+arr_el+1)/2.0
                    if (arr[arr_el]>0 and arr[arr_el+1]<0):
                        peak_indx=(arr_el+arr_el+1)/2.0
                        break
                
                a_1=Function(Q)
                a_1=interpolate(_a_trial,Q)
                arr_a1=np.flip(np.array(a_1.vector()),0)
                
                rho_1=Function(Q)
                rho_1=interpolate(_rho_trial,Q)
                arr_rho1=np.flip(np.array(rho_1.vector()),0)
                
                xstar_ind=np.argmin(abs(arr_a1-1e-3))
                xstar[n]=(xstar_ind+R_min/(-R_min+R_max)*arr_a1.size)*(-R_min+R_max)/arr_a1.size
                vstar[n]=vel(xstar[n])
                rhostar[n]=_rho_trial(xstar[n])
                
                peak_location[n]=(peak_indx+R_min/(-R_min+R_max)*arr.size)*(-R_min+R_max)/arr.size
                rho_1=Function(Q)
                rho_1=interpolate(_rho_trial,Q)
                peak_width[n]=2*(peak_indx-start_indx)*(-R_min+R_max)/arr.size
                peak_height[n]=_rho_trial(peak_location[n])-_rho_trial(start_indx*(-R_min+R_max)/arr.size)
                peak_velocity[n]=vel(peak_location[n])
                arr2=np.flip(np.array(rho_1.vector()),0)
                peak_indx=np.argmax(arr2)    
                max_location[n]=(peak_indx+R_min/(-R_min+R_max)*arr2.size)*(-R_min+R_max)/arr2.size
                arr2=np.flip(np.array(vel.vector()),0)
                indx2=np.argmax(arr2)
                max_velocity[n]=arr2[indx2]
                max_vel_loc[n]=(indx2+R_min/(-R_min+R_max)*arr2.size)*(-R_min+R_max)/arr2.size
                
                for arr_el in range(indx2-1,10,-1):
                    if (arr[arr_el]<0 and arr[arr_el+1]>0):
                        dip_location[n] = (arr_el+R_min/(-R_min+R_max)*arr.size)*(-R_min+R_max)/arr.size
                        dip_velocity[n] = (arr2[arr_el]+arr[arr_el+1])/2
                        break
                a_dip[n]=_a_trial(dip_location[n])
                rho_dip[n]=_rho_trial(dip_location[n])
                chemotaxis=project(-div(chi_0/v*grad(v)*_rho_trial),Q)
                diffusion=project(D*div(grad(rho_trial)),Q)
                growth=project(lmax*(1-_rho_trial/rhoc)*_rho_trial,Q)
            
                if (n%200)==0:
                    dum=Function(Q)
                    dum.interpolate(Expression('(tanh((4-x[0]*x[0]))+1)*.029/2',degree=3))
                    
                    plot(_rho_trial, title='Bacteria Density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'rho_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    fig=plot(chemotaxis, title='Change in density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    fig=plot(diffusion, title='Change in density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    fig=plot(growth, title='Change in density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    fig=plot(chemotaxis+diffusion+growth, title='Change in density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'drho_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.show()
                    plt.close()                    
                    
                    plot(_rho_trial, title='Bacteria Density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.yscale("log")
                    plt.ylim(0.01,np.min([rho_c,10]))
                    plt.grid()
                    plt.savefig(folder+'log_rho_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    plot(rho_g, title='Bacteria Density Gradient at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'rhog_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    plot(_rho_trial-dum, title='Difference in Bacteria Density at '+'{0:4.2f}'.format(t/60.0)+' minutes')
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'rho_diff_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    plot(_a_trial, title='Attractant Density at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'a_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    plot(_a_trial, title='Attractant Density at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.yscale("log")
                    plt.ylim(0.001,a_0)
                    plt.grid()
                    plt.savefig(folder+'log_a_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    fig=plot(_a_trial, title='Attractant/Bacteria Density at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    fig=plot(_rho_trial, title='Attractant/Bacteria Density at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.yscale("log")
                    plt.ylim(0.001,np.min([np.max([a_0,rho_c]),10]))
                    plt.grid()
                    plt.show()
                    plt.savefig(folder+'log_a_rho_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    fig=plot(_a_trial, title='Attractant/Bacteria Density/Velocity at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    fig=plot(_rho_trial, title='Attractant/Bacteria Density/Velocity at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    fig=plot(3600*vel, title='Attractant/Bacteria Density/Velocity at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.ylim(0.001,np.min([np.max([a_0,rho_c]),10]))
                    plt.yscale("log")
                    plt.grid()
                    plt.show()
                    plt.savefig(folder+'log_a_vel_rho_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                    
                    plot(vel, title='Velocity profile at '+'{0:4.2f}'.format(t/60.0)+' minutes') 
                    plt.xlabel('Distance from center of plate (in mm)')
                    plt.grid()
                    plt.savefig(folder+'vel_'+'{0:4.2f}'.format(t)+'.png', dpi=dpi)
                    plt.close()
                      
                    plt.plot(times,peak_location)
                    plt.title('Location of Peak over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1peak_loc.png', dpi=dpi)
                    plt.close()
                    
                    plt.plot(times,peak_height)
                    plt.title('Height of Peak over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1peak_hgt.png', dpi=dpi)
                    plt.close()
                    
                    plt.plot(times,peak_width)
                    plt.title('Width of Peak over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1peak_width.png', dpi=dpi)
                    plt.close()
                    
                    plt.plot(times,np.maximum(0,np.gradient(peak_location)))
                    plt.title('Expansion Speed over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1exp_speed.png', dpi=dpi)
                    plt.close()
                    
                    plt.plot(times,peak_velocity)
                    plt.title('Peak Velocity over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1peak_vel.png', dpi=dpi)
                    plt.close()
                    
                    plt.plot(times,max_velocity)
                    plt.title('Maximum Velocity over time')
                    plt.xlabel('Time (in minutes)')
                    plt.grid()
                    plt.savefig(folder+'1max_vel.png', dpi=dpi)
                    plt.close()
                    
                    scipy.io.savemat(folder+'peak_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(x=times, peak_location=peak_location))
                    scipy.io.savemat(folder+'peak_vel_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(peak_velocity=peak_velocity))
                    scipy.io.savemat(folder+'max_vel_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(max_velocity=max_velocity))
                    scipy.io.savemat(folder+'peak_width_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(peak_width=peak_width))
                    scipy.io.savemat(folder+'peak_height_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(peak_height=peak_height))
                    scipy.io.savemat(folder+'max_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(max_location=max_location))
                    scipy.io.savemat(folder+'max_vel_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(max_vel_loc=max_vel_loc))
                    scipy.io.savemat(folder+'vstar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(vstar=vstar))
                    scipy.io.savemat(folder+'xstar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(xstar=xstar))
                    scipy.io.savemat(folder+'rhostar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(rhostar=rhostar))
                    scipy.io.savemat(folder+'dip_location_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(dip_location=dip_location))
                    scipy.io.savemat(folder+'dip_velocity_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(dip_velocity=dip_velocity))
                    scipy.io.savemat(folder+'attractant_dip_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(a_dip=a_dip))
                    scipy.io.savemat(folder+'density_dip_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'.mat', dict(rho_dip=rho_dip))
                    scipy.io.savemat(folder+'velocity_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(vel_profile=arr2))
                    scipy.io.savemat(folder+'attractant_profile_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(a_profile=arr_a1))
                    scipy.io.savemat(folder+'density_profile_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(rho_profile=arr_rho1))
                    scipy.io.savemat(folder+'growth_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(growth=growth))
                    scipy.io.savemat(folder+'chemotaxis_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(chemotaxis=chemotaxis))
                    scipy.io.savemat(folder+'diffusion_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(diffusion=diffusion))
                        
                if max_vel_loc[n]-R_min>0.5*(R_max-R_min):#t>5.0*60*60:
                    #break
                    #delta=dt*1.2*(scaler*1e-6-50*1e-6)*np.sqrt(6.9/(10*60*60)*0.1/(scaler*1e-6*1e-3))
                    delta=max_vel_loc[n]-max_vel_loc[n-1]
                    print(delta)
                    R_max=R_max+delta
                    print(R_max)
                    R_min=R_min+delta#0.169/60*dt
                    print(R_min)
                    mesh = IntervalMesh(resln,R_min,R_max)
                    Q = FunctionSpace(mesh, 'P', fn_poly)
                    P1 = FiniteElement('P', interval, fn_poly)
                    element = MixedElement([P1, P1])
                    V = FunctionSpace(mesh, element)
                    rho_test, a_test = TestFunctions(V)
                    
                    # Define functions for velocity and concentrations
                    u_n = Function(V)

                u_n.assign(u)
                rho_iter, a_iter = split(u_n)
                w = Function(Q)
                u = Function(V)
                rho_trial, a_trial = split(u)
                rho_test, a_test = TestFunctions(V)
                f_1 = Constant(0.0)
                
                lam=lmax*(1-rho_trial/rhoc)
                mu=alpha*((a_trial**ahill)/((a_trial**ahill)+(Kma**ahill)))
                v=(1+a_trial/am)
                w=(chi_0*grad(a_trial)/(a_trial+am))
                F = (((rho_trial - rho_iter) / k)*rho_test*dx + div(w*rho_trial)*rho_test*dx 
                + D*dot(grad(rho_trial), grad(rho_test))*dx - lam*rho_trial*rho_test*dx  
                + ((a_trial - a_iter) / k)*a_test*dx
                + Da*dot(grad(a_trial), grad(a_test))*dx 
                + mu*rho_trial*a_test*dx
                - f_1*rho_test*dx - f_1*a_test*dx)
                
            scipy.io.savemat(upfolder+'peak_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(x=times, peak_location=peak_location))
            scipy.io.savemat(upfolder+'peak_vel_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(peak_velocity=peak_velocity))
            scipy.io.savemat(upfolder+'max_vel_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(max_velocity=max_velocity))
            scipy.io.savemat(upfolder+'peak_width_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(peak_width=peak_width))
            scipy.io.savemat(upfolder+'peak_height_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(peak_height=peak_height))
            scipy.io.savemat(upfolder+'max_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(max_location=max_location))
            scipy.io.savemat(upfolder+'max_vel_loc_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(max_vel_loc=max_vel_loc))
            scipy.io.savemat(upfolder+'vstar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(vstar=vstar))
            scipy.io.savemat(upfolder+'xstar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(xstar=xstar))
            scipy.io.savemat(upfolder+'rhostar_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(rhostar=rhostar))
            scipy.io.savemat(upfolder+'dip_location_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(dip_location=dip_location))
            scipy.io.savemat(upfolder+'dip_velocity_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(dip_velocity=dip_velocity))
            scipy.io.savemat(upfolder+'attractant_dip_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(a_dip=a_dip))
            scipy.io.savemat(upfolder+'density_dip_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(rho_dip=rho_dip))
            scipy.io.savemat(upfolder+'growth_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(growth=growth))
            scipy.io.savemat(upfolder+'chemotaxis_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(chemotaxis=chemotaxis))
            scipy.io.savemat(upfolder+'diffusion_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_'+current_time +'.mat', dict(diffusion=diffusion))
            scipy.io.savemat(folder+'velocity_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(vel_profile=arr2))            
            scipy.io.savemat(folder+'attractant_profile_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(a_profile=arr_a1))            
            scipy.io.savemat(folder+'density_profile_'+str(Da).zfill(3)+'_'+str(scaler).zfill(4)+'_{0:4.2f}'.format(t)+'.mat', dict(rho_profile=arr_rho1))
                    