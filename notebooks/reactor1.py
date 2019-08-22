def setup_initial_conditions(params):
    
    # setup the steady state for the delayed-neutron precursors
    
    n_species = len(params['species_decay'])
    
    assert len(params['species_rel_yield']) == n_species
    
    import numpy as np
    c_vec_0 = np.zeros(n_species,dtype=np.float64) # initialize conentration vector

    species_decay = params['species_decay'] # retrieve list of decay constants
    lambda_vec    = np.array(species_decay) # create a numpy vector
    
    beta = params['beta']
    species_rel_yield = params['species_rel_yield']
    beta_vec = np.array(species_rel_yield) * beta  # create the beta_i's vector

    gen_time = params['gen_time'] # retrieve neutron generation time

    n_ss = params['n_ss']
    c_vec_ss = beta_vec/lambda_vec/gen_time * n_ss # compute the steady state precursors number density
    
    params['c_vec_ss'] = c_vec_ss
    
    # setup initial condition for variables
    params['n_0']     = n_ss
    params['c_vec_0'] = c_vec_ss
    params['rho_0']   = params['reactivity']
    
    params['temp_f_0'] = 300
    params['temp_c_0'] = params['temp_0']
    params['pressure_0'] = 1.013 # bar
    
    return
def derivativeHelper(temperature):
    '''
    function to help with the implementation of alpha_tn_function since the scipy derivative function only allows
    one variable to be input to the function it works on.
    '''
    
    import iapws.iapws97 as steam
    
    pressure = pressureCalc(temperature)
    rho = 1/steam._Region1(temperature, pressure)['v']
    return(rho)

def alpha_tn_func(temp, params):
  import math
  import scipy.misc as diff
  import scipy.constants as sc
  import iapws.iapws97 as steam
    
  pressure = pressureCalc(temp)
    
  d_rho = diff.derivative(derivativeHelper, temp) # dRho/dTm
  rho = 1 / steam._Region1(temp, pressure)['v'] # mass density, kg/m3
    
  Nm = ((rho * sc.kilo)/params['mod molar mass']) * sc.N_A * (sc.centi)**3 # number density of the moderator
  d_Nm =  ((d_rho * sc.kilo)/params['mod molar mass']) * sc.N_A * (sc.centi)**3 #dNm/dTm
    
  mod_macro_a = params['mod micro a'] * Nm # macroscopic absorption cross section of the moderator
  mod_macro_s = params['mod micro s'] * Nm # macroscopic scattering cross section of the moderator
    
  F = params['fuel macro a']/(params['fuel macro a'] + mod_macro_a) # thermal utilization, F
  #dF/dTm
  d_F = -1*(params['fuel macro a'] * params['mod micro a'] * d_Nm)/(params['fuel macro a'] + mod_macro_a)**2
    
  # Resonance escape integral, P
  P = math.exp((-1 * params['n fuel'] * (params['fuel_volume']) * params['I'] * sc.zepto * sc.milli)/(mod_macro_s * params['coolant_volume']))
  #dP/dTm
  d_P = P * (-1 * params['n fuel'] * params['fuel_volume'] * sc.centi**3 * params['mod micro s'] * d_Nm)/(mod_macro_s * params['coolant_volume'] * sc.centi**3)**2
    
  Eth = 0.0862 * temp # convert temperature to energy in MeV
  E1 = mod_macro_s/math.log(params['E0']/Eth) # neutron thermalization macroscopic cross section
  Df = 1/(3 * mod_macro_s * (1 - params['mod mu0'])) # neutron diffusion coefficient
  tau = Df/E1 # fermi age, tau
    #dTau/dTm
  d_tau = (((0.0862 * (Eth/params['E0'])) * 3 * Nm) - math.log(params['E0']/Eth) * (params['mod micro s'] * d_Nm))/((3 * Nm)**2 * (1 - params['mod mu0']))
    
  L = math.sqrt(1/(3 * mod_macro_s * mod_macro_a * (1 - params['mod mu0']))) # diffusion length L
    # dL/dTm
  d_L = 1/(2 * math.sqrt((-2 * d_Nm)/(3 * params['mod micro s'] * params['mod micro a'] * Nm**3 * (1 - params['mod mu0']))))
    
  # left term of the numerator of the moderator temperature feedback coefficient, alpha
  left_1st_term = d_tau * (params['buckling']**2 + L**2 * params['buckling']**4) #holding L as constant
  left_2nd_term = d_L * (2 * L * params['buckling']**2 + 2 * L * tau * params['buckling']**4) # holding tau as constant
  left_term = (P * F) * (left_1st_term + left_2nd_term) # combining with P and F held as constant
    
  # right term of the numerator of the moderator temperature feedback coefficient, alpha
  right_1st_term = (-1) * (1 + ((tau + L**2) * params['buckling']**2) + tau * L**2 * params['buckling']**4) # num as const
  right_2nd_term = F * d_P # holding thermal utilization as constant
  right_3rd_term = P * d_F # holding resonance escpae as constant
  right_term = right_1st_term * (right_2nd_term + right_3rd_term) # combining all three terms together
    
  # numerator and denominator
  numerator = left_term + right_term
  denominator = params['eta'] * params['epsilon'] * (F * P)**2
  alpha_tn = numerator/denominator
    
  alpha_tn = alpha_tn * (sc.milli * sc.zepto * 0.000000001) # adjust for barns
  return alpha_tn

def rho_func( t, n_dens, temp, params ):
  '''
  Reactivity function.  
    
  Parameters
  ----------
  t: float, required
      Time.
  temp_f: float, required
      Temperature at time t.
  params: dict, required
      Dictionary of quantities. It must have a `'rho_0'` key/value pair.
    
  Returns
  -------
  rho_t: float
      Value of reactivity.

  Examples
  --------
  '''
    
  rho_0  = params['rho_0']
  temp_ref = params['temp_c_ss_operation']
  n_dens_ss_operation = params['n_dens_ss_operation']
  alpha_n = params['alpha_n']
    
  alpha_tn = alpha_tn_func(temp, params)
    
  if t > params['shutdown time']:
      rho_0 = 0
    
  rho_t = rho_0 + alpha_n * n_dens + alpha_tn * (temp - temp_ref)
  return rho_t
  
def q_source( t, params ):
  '''
  Neutron source delta function.  
    
  Parameters
  ----------
  t: float, required
      Time.
  params: dict, required
      Dictionary of quantities. It must have a `'q_0'` key/value pair.
    
  Returns
  -------
  q: float
      Value of source.

  Examples
  --------
  '''
    
  q = 0
  q_0 = params['q_0']
    
  if t <= 1e-5: # small time value
      q = q_0
  else:
      q = 0.0
        
  return q
  
def sigma_fis_func( temp, params ):
  '''
  Place holder for implementation
  '''
  import math
  sigma_f = params['sigma_f_o'] * math.sqrt(298/temp) * math.sqrt(math.pi) * 0.5
    
  return sigma_f
  
def nuclear_pwr_dens_func( time, temp, n_dens, params ):
  '''
  Place holder for implementation
  '''
  rxn_heat = params['fis_energy'] # get fission reaction energy J per reaction
    
  sigma_f = sigma_fis_func( temp, params ) # m2
    
  fis_nuclide_num_dens = params['fis_nuclide_num_dens_fake'] #  #/m3
    
  Sigma_fis = sigma_f * fis_nuclide_num_dens # macroscopic cross section
    
  v_o = params['thermal_neutron_velo'] # m/s
    
  neutron_flux = n_dens * 9.08E15 * v_o
    
   #reaction rate density
  rxn_rate_dens = Sigma_fis * neutron_flux
    
  # nuclear power source
  q3prime = - rxn_heat * rxn_rate_dens # exothermic reaction W/m3
  #q3prime = - n_dens * 3323E6
  #print("q3prime")
  #print(q3prime)
    
  return q3prime
  
def heat_sink_rate( time, temp_f, temp_c, params):
    
  ht_coeff = params['ht_coeff']
    
  q_f = - ht_coeff * (temp_f - temp_c)
  #print(q_f)
  return q_f
  
def pressureCalc(temperature):
  # calculates the pressure exerted by steam in the reactor at any temperature T. Returns units of mPa
  import scipy.constants as const
  n = 101111 # moles
  V = 150 # m^3
    
  pressure = n * const.R * temperature / V # PV = NRT
    
  #print(pressure)
    
  return pressure/const.mega
  
def f_vec(time, u_vec, params):
  import numpy as np
  assert np.all(u_vec >= 0.0)
    
  n_dens = u_vec[0] # get neutron dens

  c_vec = u_vec[1:-2] # get delayed neutron emitter concentration
    
  temp_f = u_vec[-2] # get temperature of fuel
    
  temp_c = u_vec[-1] # get temperature of coolant
    
  # initialize f_vec to zero
  species_decay = params['species_decay']
  lambda_vec = np.array(species_decay)
  n_species  = len(lambda_vec)
    
  f_tmp = np.zeros(1+n_species+2,dtype=np.float64) # vector for f_vec return
    
  #----------------
  # neutron balance
  #----------------
  rho_t    = rho_func(time, n_dens, (temp_f+temp_c)/2.0, params)
    
  beta     = params['beta']
  gen_time = params['gen_time']
       
  species_rel_yield = params['species_rel_yield']
  beta_vec = np.array(species_rel_yield) * beta
    
  assert len(lambda_vec)==len(beta_vec)
    
  q_source_t = q_source(time, params)
    
  f_tmp[0] = (rho_t - beta)/gen_time * n_dens + lambda_vec @ c_vec + q_source_t
    
  #-----------------------------------
  # n species balances (implicit loop)
  #-----------------------------------
  f_tmp[1:-2] = beta_vec/gen_time * n_dens - lambda_vec * c_vec
    
  #--------------------
  # fuel energy balanc
  #--------------------
  rho_f    = params['fuel_dens']
  cp_f     = params['cp_fuel']
  vol_fuel = params['fuel_volume']
    
  pwr_dens = nuclear_pwr_dens_func( time, (temp_f+temp_c)/2, n_dens, params )
    
  heat_sink = heat_sink_rate( time, temp_f, temp_c, params )
  #assert heat_sink <= 0.0,'heat_sink = %r'%heat_sink
    
  f_tmp[-2] =  -1/rho_f/cp_f * ( pwr_dens - heat_sink/vol_fuel )
    
  #-----------------------
  # coolant energy balance
  #-----------------------
  rho_c    = params['coolant_dens']
  cp_c     = params['cp_coolant']
  vol_cool = params['coolant_volume']
    
  # subcooled liquid
  turbine_out = turbine(time, temp_c, params)[0] #run the turbine, take the runoff and pass to condenser
  condenser_out = condenser(time, turbine_out, temp_c, params)[0] #run the condenser, pass runoff to the pump
  pump_out = pump(time, condenser_out, temp_c, params) #run the pump, runoff returns to reactor as temp_in
  #print("time is ", time, "and inlet temperature is", temp_in, "\n")
    
  temp_in = pump_out
    
  tau = params['tau_fake']
    
  heat_source = - heat_sink
    
  f_tmp[-1] = - 1/tau * (temp_c - temp_in) + -1./rho_c/cp_c/vol_cool * heat_source
    
  # pressure calculations

  #print(time)
  #print(u_vec)
  return f_tmp
  
def run_point_reactor( f_vec, params ):
  import numpy as np
  from scipy.integrate import odeint # Load ODE solver package

  import numpy as np
  time_final = params['time_final']
  n_time_stamps = params['n_time_stamps']
  time_stamps = np.linspace(0.0, time_final, num=n_time_stamps) # create the time stamps for solution values
  params['time_stamps'] = time_stamps
    
  max_n_steps_per_time_step = 1000 # max number of nonlinear algebraic solver iterations per time step

  n_0     = params['n_0']
  c_vec_0 = params['c_vec_0']
    
  temp_f_0 = params['temp_f_0']
  temp_c_0 = params['temp_c_0']
  pressure_0 = params['pressure_0']
       
  # m-equation point reactor model
  n_species = len(c_vec_0)
  u_vec_0 = np.zeros(1+n_species+2,dtype=np.float64)
    
  u_vec_0[0]    = n_0
  u_vec_0[1:-2] = c_vec_0
  u_vec_0[-2]   = temp_f_0
  u_vec_0[-1]   = temp_c_0
  #u_vec_0[-1]   = pressure_0
            
  (u_vec_history, info_dict) = odeint( f_vec, u_vec_0, time_stamps,
                                       args=( params, ),
                                       rtol=1e-4, atol=1e-8, mxstep=max_n_steps_per_time_step,
                                       full_output=1, tfirst=True )
  #print(n_dens)
  #print(u_vec)
  #print(time_stamps)
  #print

  assert info_dict['message']=='Integration successful.',\
                   'Fatal: scipy.integrate.odeint failed %r'%info_dict['message']
    
  return u_vec_history
  
def plot_results( u_vec_history, params, normalize=True, semi_log=False, markers=False, precursors=True ):
  import numpy as np
    
  time_stamps = params['time_stamps']/3600
  tau = params['tau_fake']
        
  import matplotlib.pyplot as plt
    
  fig, ax1 = plt.subplots(1, figsize=(14, 6))

  if precursors == True:
        
      ax2 = ax1.twinx() # duplicate x axes to plot n and c_i's in different y axes
    
      color_ids = np.linspace(0,1,u_vec_history[:,1:-2].shape[1])
    
      for (j,color_id) in zip( range(u_vec_history[:,1:-2].shape[1]), color_ids ):
          color=plt.cm.nipy_spectral(color_id)
        
          if normalize == True:
              ax2.plot( time_stamps,u_vec_history[:,j+1]/params['c_vec_0'][j],'-.',color=color,label=r'$c_%i$'%(j+1) )
              ax2.set_ylabel(r'$c_i/c_{i_0}$',fontsize=16,color='black')
          else:
              ax2.plot( time_stamps,u_vec_history[:,j+1],'-.',color=color,label=r'$c_%i$'%(j+1) )
              ax2.set_ylabel(r'$c_i$',fontsize=16,color='black')
        
      ax2.tick_params(axis='y', labelcolor='black', labelsize=14)
      ax2.legend(loc='lower right',fontsize=12)
      if semi_log == True:
          ax2.set_yscale('log') # uncomment to plot y in log scale
      #ax2.grid(True)

  if markers == True:
      if normalize == True:
          ax1.plot( time_stamps,u_vec_history[:,0]/params['n_0'],'-',marker='+',color='red',label=r'$n/n_0$' )
          ax1.set_ylabel(r'$n$',fontsize=16,color='black')
      else:
          ax1.plot( time_stamps,u_vec_history[:,0],'-',marker='+',color='red',label=r'$n$' )
          ax1.set_ylabel(r'$n$',fontsize=16,color='black')
  else:
      if normalize == True:
          ax1.plot(time_stamps,u_vec_history[:,0]/params['n_0'],'-',color='red',label=r'$n/n_0$' )
          ax1.set_ylabel(r'$n/n_0$',fontsize=16,color='black')
      else:
          ax1.plot(time_stamps,u_vec_history[:,0],'-',color='red',label=r'$n$' )
          ax1.set_ylabel(r'$n$',fontsize=16,color='black')

  ax1.set_xlabel(r'Time [h]',fontsize=16)
    
  ax1.tick_params(axis='y', labelcolor='black', labelsize=14)
  ax1.tick_params(axis='x', labelsize=14)
  ax1.legend(loc='best',fontsize=12)
  if semi_log == True:
      ax1.set_yscale('log') # uncomment to plot y in log scale
  ax1.grid(True)

  plt.title(r'Point-Reactor Model: $\rho/\beta=$'
            +str(params['reactivity']/params['beta'])
            +r'; $q_0=$'+str(round(params['q_0'],2)),
            fontsize=18)
  fig.tight_layout()  # otherwise the right y-label is slightly clipped
  plt.show() 

  print('')
  return
  
def peek(time,data, head=500, tail=100):  

  import pandas as pd
    
  pd.options.display.float_format = '{:.2e}'.format
    
  layout = {'time':time[:head]}
    
  layout['n'] = data[:head,0]
    
  for j in range(1,data[:,1:-2].shape[1]+1):
      layout['c_%i'%j] = data[:head,j]
        
  layout['temp_f'] = data[:head,-2]
  layout['temp_c'] = data[:head,-1]
        
  results = pd.DataFrame(layout)
  print(round(results,2))
  print('')
    
  #layout = {'time':time[-tail:]}
    
  #layout['n'] = data[-tail:,j]
  #for j in range(1,data[:,1:-2].shape[1]+1):
  # layout['c_%i'%j] = data[-tail:,j]
    
  # layout['temp_f'] = data[-tail:,-2]
 # layout['temp_c'] = data[-tail:,-1]
               
 # results = pd.DataFrame(layout)
 # print(round(results,2))
 # print('')
  return
  
def turbine (time, temp_in, params):
  import iapws.iapws97 as steam_table
  #expand the entering steam from whatever temperature and pressure it enters at to 0.035 kpa, with 80% efficiency.
  #pressure of steam when it enters the turbine equals the current reactor operating pressure
  pressure = pressureCalc(temp_in) 
        
  #print("Pressure is " + str(pressure) + " and temperature is " + str(temp_in))
        
      #enthalpy of the steam coming out of the reactor; for now, assume 100% quality or superheated steam
  h_steam = steam_table._Region1(temp_in, pressure)["h"] 
  h_liquid = steam_table._Region4(0.00075, 0)["h"]#enthalpy of the ideal liquid runoff
        
      #print(pressure, temp_in)
        
  w_isentropic = h_steam - h_liquid # isentropic work
  w_real = params['turbine efficiency'] * w_isentropic # actual work
  h_real = h_steam - w_real # h_end = h_in - w
  w_real = w_real * params['steam flowrate'] #multiply work done/kg steam by amount of steam to get total work
            
  t_runoff = steam_table._Backward1_T_Ph(0.00075, h_real) # this goes to the condenser
  return (t_runoff, w_real)
  
def condenser(time, temp_in, temp_c, params):
  #compress the liquid to 101.3 kpa and store the work done by the condenser
  import iapws.iapws97 as steam_table
  pressure = pressureCalc(temp_c) # current reactor operating pressure that the runoff must be compressed to
    
  subcooling = params['% subcooling']
    
  #condenser functionality is partially compromised, leading to a lower degree of exit subcooling
  if time > params['malfunction start'] and time < params['malfunction end']:
      subcooling = params['malfunction subcooling']
    
  # condenser fails completely, and as a result coolant flow to and from the core is compromised
  #elif time > params['breakage start'] and time < params['breakage end']:
      #params['steam flowrate'] = params['breakage steam flowrate']
      #params['coolant_volume'] = params['breakage coolant_volume']
      #params['breakage reached'] = True
        
  #if time > params['shutdown time'] and params['shutdown temp reached'] == False:
      #params['% subcooling'] = 0.99 * params['% subcooling']
    
  # during shutdown, the degree of subcooling is gradually increased
        
 # if params['breakage reached'] == True and time > params['breakage end']:
  #    params['steam flowrate'] = params['normal steam flowrate']
  #    params['coolant_volume'] = params['normal coolant_volume']
    #  params['breakage reached'] = False
    
  h_saturated = steam_table._Region4(pressure, 0)["h"] # enthalpy of water at 1 atm
  h_subcooled = subcooling * h_saturated #enthalpy of the desired liquid runoff from the condenser
        
  h_turbine = steam_table._Region1(temp_in, 0.00075)["h"] #enthalpy of the liquid leaving the turbine at 0.035 kPa
        
  work_done = h_subcooled - h_turbine #work done by the condenser per kg steam effluent
  work_done = work_done * params['steam flowrate'] / params['condenser efficiency']
        
  t_runoff = steam_table._Backward1_T_Ph(pressure, h_subcooled) #runoff temperature
    
  #ensure that shutdown temp hasn't been breached
  #if t_runoff <= 310:
      #t_runoff = 310
      #params['% subcooling'] = params['% subcooling']/0.99
      #params['shutdown temp reached'] = True
  return (t_runoff, work_done)
  
def pump(time, temp_effluent, temp_c, params):
  #placeholder for a more accurate pumping function that includes a slight temperature increase and work done
  return temp_effluent
  
def quantities1(u_vec_history, params, time_stamps):
  # used to graph q''' and heat removed
  import pandas as pd
  import scipy.constants
  data = dict()
    
  q3_list = list() #temp storage for q3 prime
  removed_heat_list = list() #manipulated q3 prime
    
  for (time, n_dens, temp_f, temp_c) in zip(time_stamps, u_vec_history[:,0], u_vec_history[:,-2], u_vec_history[:,-1]):
       
      q3prime = abs(nuclear_pwr_dens_func(time, temp_f, n_dens, params)) # calculate q3prime at this point in time; watts/m3
      q3_list.append(q3prime/scipy.constants.giga) # convert watts/m3 to gWatts/m3
        
      heat_removed = -1 * heat_sink_rate("false", temp_f, temp_c, params) # calculate the heat removed at this point in time
      removed_heat_list.append(heat_removed/scipy.constants.giga) # convert gwatts/m3 to kWatts/m3
        
  #data['time [s]'] = time_stamps
  data["q''' [gW/m3]"] =  q3_list 
  data["heat removed [gW/m3]"] = removed_heat_list
    
  quantities = pd.DataFrame( data )
  return(quantities)
  
def quantities2(u_vec_history, params, time_stamps):
  # used to graph turbine, condenser and pump work
  data = dict()
  import pandas as pd
  import scipy.constants
    
  twork = list() #turbine work
  cwork = list() #condenser work
  nwork = list() #net total usable work the reactor and BOP produce
  #pwork = list() # pump work
    
  for (time, n_dens, temp_f, temp_c) in zip(time_stamps, u_vec_history[:,0], u_vec_history[:,-2], u_vec_history[:,-1]):
        
      turb = turbine(time, temp_c, params)
      turbwork = turb[1]/scipy.constants.mega
      twork.append(turbwork)
        
      cond = condenser(time, turb[0], temp_c, params)
      condwork = cond[1]/scipy.constants.mega
      cwork.append(condwork)
        
      net_work = turbwork - abs(condwork) # Wnet = Ws - |Wcond|
      nwork.append(net_work)
        
      #pwork.append(pump_work[time])
        
  data["turbine work [mW]"] = twork
  data["condenser work [mW]"] = cwork
  data["net work [mW]"] = nwork
    
  quantities = pd.DataFrame( data )
  return(quantities)
  
def tmp(u_vec_history, params):
  time_stamps = params['time_stamps']
  tau = params['tau_fake']
  import matplotlib.pyplot as plt    
  fig, ax1 = plt.subplots(1, figsize=(16, 6))
  ax1.plot(time_stamps/tau,u_vec_history[:,-2],'b-',label='$T_f=$ ' )

  ax1.set_xlabel(r'Time [s] ($\tau=$%4.1f s)'%tau,fontsize=16)
  ax1.set_ylabel(r'$T_f$ [K]',fontsize=16,color='blue')
  ax1.tick_params(axis='y', labelcolor='blue', labelsize=14)
  ax1.tick_params(axis='x', labelsize=14)
  ax1.legend(loc='best',fontsize=12)
  ax1.grid(True)

  ax2 = ax1.twinx() 
  ax2.plot(time_stamps/tau,u_vec_history[:,-1],'g-.',label='$T_c=$ ' )
  ax2.set_ylabel(r'$T_c$ [K]',fontsize=16,color='green')
  ax2.tick_params(axis='y', labelcolor='green', labelsize=14)
  ax2.legend(loc='best',fontsize=12)
  #ax2.grid(True)

  plt.title('CSTR w/ Cooling Coil (exothermic rxn)',fontsize=20)
  fig.tight_layout()  # otherwise the right y-label is slightly clipped
  plt.show()
  print('')
  return