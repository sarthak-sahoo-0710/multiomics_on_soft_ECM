#%%

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16, 'font.family': 'Arial'})

# =============================================================================
# === 1. MODEL DEFINITION (ODE FUNCTION AND EVENTS)
# =============================================================================
import numpy as np

def yap_cell_cycle_model(t, y, p, E):
    """
    MERGED MODEL (Time unit: hours)
    
    All parameters are pre-scaled to /h in the parameter dictionary.
    Both Model 1 (YAP/TAZ) and Model 2 (Cell Cycle) derivatives are calculated 
    directly in /h and used as-is.
    """
    
    # --- Unpack All 39 Species ---
    # Model 1 Species (Indices 0-26)
    FAK, pFAK, RhoA_GDP, RhoA_GTP, ROCK, ROCK_A, mDia, mDia_A, Myo, Myo_A, \
    LIMK, LIMK_A, Cofilin_p, Cofilin_NP, G_actin, F_actin, YAPTAZ_p, YAPTAZ, \
    LaminA_p, LaminA, NPC, NPC_A, YAPTAZ_nuc, Importin, mybl2, yap_mybl2, \
    SF_remod = y[:27]
    
    # Model 2 Species (Indices 27-39)
    Bm, CycB, C20m, Cdc20t, cdhm, cdht, Cdh1, Cdc20A, IEP, cdt1, EpiSil, Ecad_m, Ecad = y[27:]
    
    # =================================================================
    # --- MODULE 1: YAP/TAZ Mechanotransduction (Rates calculated in /h) ---
    # =================================================================

    basal_activation = p['k_f'] * FAK # basal FAK activation 
    stiffness_activation = p['k_sf'] * (E / (p['C'] + E)) * FAK # stiffness-dependent FAK activation
    R1 = basal_activation + stiffness_activation # total FAK activation rate
    R2 = p['k_df'] * pFAK # FAK deactivation rate
    
    rho_activation_rate = p['k_fk_rho'] * (p['gamma'] * pFAK**p['n'] + 1) # RhoA activation via FAK
    mybl2_feedback_rate = p['k_rho_mybl2'] * yap_mybl2  # RhoA activation via mybl2 feedback
    R3 = (rho_activation_rate + mybl2_feedback_rate) * RhoA_GDP # RhoA activation rate
    R4 = p['k_d_rho'] * RhoA_GTP # RhoA deactivation rate

    R5 = p['k_drock'] * ROCK_A # ROCK deactivation rate
    R6 = p['k_r_rho'] * RhoA_GTP * ROCK # ROCK activation rate
    T_ROCKA = (np.tanh(p['sc1'] * (ROCK_A - p['ROCKb'])) + 1) * ROCK_A / 2 # Thresholded active ROCK

    R7 = p['k_dmdia'] * mDia_A # mDia deactivation rate
    R8 = p['k_m_rho'] * RhoA_GTP * mDia # mDia activation rate

    R9 = p['k_mr'] * (p['epsilon'] * T_ROCKA + 1) * Myo - p['k_dmy'] * Myo_A # Myosin activation/deactivation
    R10 = p['k_lr'] * (p['tau'] * T_ROCKA + 1) * LIMK - p['k_dl'] * LIMK_A # LIMK activation/deactivation
    R11 = p['k_turn_over'] * Cofilin_p - (p['k_catcofilin'] * LIMK_A * Cofilin_NP) / (p['k_mcofilin'] + Cofilin_NP) # Cofilin phosphorylation/dephosphorylation
    R12 = p['k_ra'] * (p['alpha'] * T_ROCKA + 1) * G_actin - (p['k_dep'] + p['k_fc1'] * Cofilin_NP) * F_actin # Actin polymerization/depolymerization
    
    R13 = (p['k_CN'] + p['k_CY'] * F_actin * Myo_A) * YAPTAZ_p - p['k_NC'] * YAPTAZ # YAPTAZ phosphorylation/dephosphorylation
    E_cytosol = p['p'] * F_actin**2.6 # Cytosolic stiffness due to stress fibers
    R14 = p['k_fl'] * (E_cytosol / (p['C_Lamin'] + E_cytosol)) * LaminA_p - p['k_rl'] * LaminA # LaminA phosphorylation/dephosphorylation
    R15 = p['k_fNPC'] * LaminA * F_actin * Myo_A * NPC - p['k_rNPCA'] * NPC_A # NPC activation/deactivation
    R16 = (p['k_inb_max'] * (Importin/3) + p['k_in'] * NPC_A) * YAPTAZ - p['k_out'] * YAPTAZ_nuc # YAPTAZ nuclear import/export 
    
    production_imp = p['basal_imp_prod'] + p['k_prod_imp'] * (YAPTAZ_nuc**p['h_imp'] / (p['K_imp']**p['h_imp'] + YAPTAZ_nuc**p['h_imp'])) # importin production is modeled as a function of nuclear YAP/TAZ
    degradation_imp = p['k_deg_imp'] * Importin # importin degradation
    dImportin_dt_h = production_imp - degradation_imp # net production of importin
    
    production_mybl2 = p['basal_mybl2_prod'] + p['k_prod_mybl2'] * (YAPTAZ_nuc**p['h_mybl2'] / (p['K_mybl2']**p['h_mybl2'] + YAPTAZ_nuc**p['h_mybl2'])) # mybl2 production as a function of nuclear YAP/TAZ
    degradation_mybl2 = p['k_deg_mybl2'] * mybl2 # mybl2 degradation
    net_prod_mybl2 = production_mybl2 - degradation_mybl2 # net production of mybl2

    d_yap_mybl2_dt_h = p['comp_f_ym'] * mybl2 * YAPTAZ_nuc - p['k_diss_ym'] * yap_mybl2 # formation/dissociation of yap-mybl2 complex

    k_remod_f_eff = p['k_remod_f'] if E < p['E_sil_thresh'] else 0.0
    k_remod_d_eff = p['k_remod_d'] if E < p['E_sil_thresh'] else 0.0
    dSF_remod_dt_h = k_remod_f_eff * F_actin - k_remod_d_eff * SF_remod

    production_FAK = p['FAK_prod_rate']

    dFAK_dt_h = R2 - R1 + production_FAK*(1 - (FAK+pFAK)/(p['fak_upreg_K']))* SF_remod - p['deg_FAK']*FAK* SF_remod # net FAK dynamics with production and degradation linked to SF_remod
    dpFAK_dt_h = R1 - R2 # pFAK dynamics
    dRhoA_GDP_dt_h = -R3 + R4 # RhoA dynamics
    dRhoA_GTP_dt_h = R3 - R4
    dROCK_dt_h = R5 - R6 # ROCK dynamics
    dROCK_A_dt_h = -R5 + R6
    dmDia_dt_h = R7 - R8 # mDia dynamics
    dmDia_A_dt_h = -R7 + R8
    dMyo_dt_h = -R9 # Myosin dynamics
    dMyo_A_dt_h = R9
    dLIMK_dt_h = -R10 # LIMK dynamics
    dLIMK_A_dt_h = R10
    dCofilin_p_dt_h = -R11 # Cofilin dynamics
    dCofilin_NP_dt_h = R11
    dG_actin_dt_h = -R12 # Actin dynamics
    dF_actin_dt_h = R12
    dYAPTAZ_p_dt_h = -R13 # YAPTAZ dynamics
    dYAPTAZ_dt_h = R13 - R16
    dLaminA_p_dt_h = -R14 # LaminA dynamics
    dLaminA_dt_h = R14
    dNPC_dt_h = -R15 # NPC dynamics
    dNPC_A_dt_h = R15
    dYAPTAZ_nuc_dt_h = R16 - d_yap_mybl2_dt_h # Nuclear YAPTAZ dynamics
    dmybl2_dt_h = net_prod_mybl2 - d_yap_mybl2_dt_h # mybl2 dynamics
    dyap_mybl2_dt_h = d_yap_mybl2_dt_h # yap-mybl2 complex dynamics

    # List of Model 1 derivatives (in /h)
    derivs_model1_h = [
        dFAK_dt_h, dpFAK_dt_h, dRhoA_GDP_dt_h, dRhoA_GTP_dt_h, dROCK_dt_h, dROCK_A_dt_h,
        dmDia_dt_h, dmDia_A_dt_h, dMyo_dt_h, dMyo_A_dt_h, dLIMK_dt_h, dLIMK_A_dt_h,
        dCofilin_p_dt_h, dCofilin_NP_dt_h, dG_actin_dt_h, dF_actin_dt_h, dYAPTAZ_p_dt_h,
        dYAPTAZ_dt_h, dLaminA_p_dt_h, dLaminA_dt_h, dNPC_dt_h, dNPC_A_dt_h,
        dYAPTAZ_nuc_dt_h, dImportin_dt_h, dmybl2_dt_h, dyap_mybl2_dt_h, dSF_remod_dt_h
    ]

    # =================================================================
    # --- MODULE 2: Cell Cycle (Rates calculated in /h) ---
    # =================================================================

    # sanity checks to avoid division by zero
    cdht_minus_cdh1 = max(0, cdht - Cdh1)
    denom_cdh1 = (p['J3']) + cdht_minus_cdh1
    if denom_cdh1 == 0: denom_cdh1 = 1e-9

    cdc20t_minus_cdc20a = max(0, Cdc20t - Cdc20A)
    denom_cdc20a = (p['J7']) + cdc20t_minus_cdc20a
    if denom_cdc20a == 0: denom_cdc20a = 1e-9

    dBm_dt_h = (((p['k1m'] * p['d']) * p['GF']) / (p['kmm'] + (p['keff'] * p['GF'])) - (p['k1dm'] * p['d']) * Bm)  # Dynamics for the cyclin B mRNA
    dCycB_dt_h = (p['k1'] * Bm - p['k2a'] * CycB - p['k2b'] * CycB * Cdh1) * p['d'] # Dynamics for the cyclin B protein

    cycb_hill = (CycB ** p['n_cc']) / ((p['J5']) ** p['n_cc'] + (CycB ** p['n_cc'])) # Hill function for cyclin B effect on cdc20 mRNA production
    mybl2_activation = (yap_mybl2**p['h_mybl2_cdc']) / (p['K_mybl2_cdc']**p['h_mybl2_cdc'] + yap_mybl2**p['h_mybl2_cdc']) # Hill function for mybl2 effect on cdc20 mRNA production
    
    basal_prod_c20m = (p['k5am']) 
    cycb_prod_c20m = (p['k5bm'] * cycb_hill) / (p['k5cm'] + (p['GF'] * p['j5c']))
    total_prod_c20m = (basal_prod_c20m + cycb_prod_c20m) * mybl2_activation # Total production rate of cdc20 mRNA
    dC20m_dt_h = (total_prod_c20m - p['k5dm'] * C20m) * p['d'] # Dynamics for the cdc20 mRNA

    dCdc20t_dt_h = (p['k5a'] * C20m - p['k6'] * Cdc20t) * p['d'] # Dynamics for the total cdc20 protein

    # Base calculations for silencing and surge
    silencing_factor = 1.0 / (1.0 + (SF_remod / p['k_sil_strength'])**p['hills_coeff_silencing'])
    ecad_normalized = max(0.0, (Ecad / p['Ecad_ref_scale']) - 1.0)
    fzr1_hill_surge = (ecad_normalized**p['hills_coeff_fzr1_surge']) / (ecad_normalized**p['hills_coeff_fzr1_surge'] + (p['Fzr1_Ecad_thres'])**p['hills_coeff_fzr1_surge'])

    # Apply model specific overrides based on 'model_type' parameter (default to 3 for base model)
    model_type = p.get('model_type', 3)
    if model_type == 1:
        fzr1_hill_surge = 0.0
    elif model_type == 2:
        silencing_factor = 1.0

    k3m_effective = p['k3m'] * (1.0 + p['upreg_k3m_fzr1'] * fzr1_hill_surge) if E < p['E_sil_thresh'] else p['k3m']

    dcdhm_dt_h = (k3m_effective - p['k3dm'] * cdhm) * p['d'] # dynamics for cdh1 mRNA
    dcdht_dt_h = (p['k3a'] * cdhm - p['k3dt'] * cdht) * p['d'] # dynamics for cdh1 total protein

    ke_effective = p['k3'] 
    term1_cdh1 = ((ke_effective) + (p['k3b'] * Cdc20A)) * cdht_minus_cdh1 / denom_cdh1 # activation term for cdh1
    term2_cdh1 = p['k4'] * CycB * Cdh1 / ((p['J4']) + Cdh1) # inactivation term for cdh1
    dCdh1_dt_h = (term1_cdh1 - term2_cdh1 - p['k3dt'] * Cdh1) * p['d'] # dynamics for cdh1 active protein

    term1_cdc20a = p['k7'] * IEP * cdc20t_minus_cdc20a / denom_cdc20a
    term2_cdc20a = p['k8'] * p['Mad'] * Cdc20A / ((p['J8']) + Cdc20A)
    dCdc20A_dt_h = (term1_cdc20a - term2_cdc20a - p['k6'] * Cdc20A) * p['d'] # dynamics for active cdc20 protein

    dIEP_dt_h = (p['k9'] * (CycB) * (1 - IEP) - p['k10'] * IEP) * p['d'] # dynamics for IEP
    dcdt1_dt_h = (p['k11'] - p['k12'] * CycB * cdt1 - p['k13'] * cdt1) * p['d'] # dynamics for cdt1

    dEpiSil_dt_h = p['k_epi_on'] * ((SF_remod**3)/(p['threshold_remod_to_DNMT_act']**3 + SF_remod**3)) - p['Epi_deg'] * EpiSil  #this equation and associated parameters have not been used in the model

    k_surge_ecad_prod = (p['k_ecad_prod'] + ((p['Ecad_upreg']*p['k_ecad_prod'])*silencing_factor)) if E < p['E_sil_thresh'] else p['k_ecad_prod']
    dEcad_m_dt_h = (k_surge_ecad_prod  - p['k_ecad_deg_m'] * Ecad_m) * p['d']
    dEcad_dt_h = (p['k_ecad_trans'] * Ecad_m - p['k_ecad_deg_p'] * Ecad) * p['d']

    derivs_model2_h = [
        dBm_dt_h, dCycB_dt_h, dC20m_dt_h, dCdc20t_dt_h, dcdhm_dt_h, dcdht_dt_h,
        dCdh1_dt_h, dCdc20A_dt_h, dIEP_dt_h, dcdt1_dt_h, dEpiSil_dt_h, dEcad_m_dt_h, dEcad_dt_h
    ]
    
    return derivs_model1_h + derivs_model2_h

def division_event(t, y, p, E):
    return y[28] - 0.1 # trigger when cycb levels fall

division_event.terminal = True  # Stop integration at event
division_event.direction = -1  # Trigger only when decreasing

# =============================================================================
# === 2. SETUP FUNCTIONS (PARAMETERS AND INITIAL CONDITIONS)
# =============================================================================

def get_default_parameters():
    # --- Define Model 1 Parameters (units are pre-scaled to /h) ---
    params_model1_hr = {
        'k_f': 0.015*3600, 'k_sf': 0.379*3600, 'k_df': 0.035*3600, 'C': 3.25, 'k_fk_rho': 0.0168*3600,
        'gamma': 77.56, 'n': 5, 'k_d_rho': 0.625*3600, 'k_r_rho': 0.648*3600, 'k_drock': 0.8*3600,
        'k_m_rho': 0.002*3600, 'k_dmdia': 0.005*3600, 'k_mr': 0.03*3600, 'epsilon': 36, 'k_dmy': 0.067*3600,
        'sc1': 20, 'ROCKb': 0.3, 'k_lr': 0.07*3600, 'tau': 55.49, 'k_dl': 2*3600, 'k_turn_over': 0.04*3600,
        'k_catcofilin': 0.34*3600, 'k_mcofilin': 4*3600, 'k_ra': 0.4*3600, 'alpha': 50, 'k_dep': 3.5*3600,
        'k_fc1': 4*3600, 'k_CN': 0.56*3600, 'k_CY': 7.6e-4*3600, 'k_NC': 0.14*3600, 'k_in': 10*3600, 'k_out': 1*3600,
        'k_fl': 0.46*3600, 'C_Lamin': 100, 'k_rl': 0.001*3600, 'k_fNPC': 2.8e-7*3600, 'k_rNPCA': 8.7*3600, 'p': 9e-6,
        'k_inb_max': 1.0*3600, 
        'k_prod_imp': 0.39*3600, 'k_deg_imp': 1*3600, 'K_imp': 0.4, 'h_imp': 4, 'basal_imp_prod': 2.61*3600,
        'k_prod_mybl2': 0.2*3600, 'k_deg_mybl2': 0.1*3600, 'K_mybl2': 0.4, 'h_mybl2': 4, 'basal_mybl2_prod': 0.01*3600,
        'comp_f_ym': 0.3*3600, 'k_diss_ym': 0.53*3600,
        'k_rho_mybl2': 0.0015*3600,
        'k_remod_f': 1e-7*3600, 'k_remod_d': 1e-5*3600, 'FAK_prod_rate': 0.35*3600, 'deg_FAK': 0.15*3600, 'fak_upreg_K': 1.25, 'threshold_remod_to_DNMT_act': 0.5
    }
    
    # --- Define Model 2 Parameters (units are /h) ---
    params_model2_hr = {
        'GF': 2.0, 'kmm': 0.20, 'keff': 1.0, 'sf': 1.0, 'd': 2.8,
        'k1m': 0.0037, 'k1dm': 0.058, 'k1': 0.4, 'k2a': 0.04, 'k2b': 2.0,
        'k3': 1.28, 'k3a': 1.0, 'k3b': 8.0, 'k4': 40.0, 'J3': 0.04, 'J4': 0.04,
        'k3m': 0.5, 'k3dm': 0.5, 'k3dt': 1.0,
        'k5am': 0.005, 'k5bm': 0.2, 'k5cm': 1.0, 'j5c': 0.02, 'k5dm': 1.386,
        'J5': 0.3, 'n_cc': 4, 'k5a': 1.0, 'k6': 0.05, 
        'k7': 1.4, 'k8': 0.5, 'J7': 0.001, 'J8': 0.001, 'Mad': 1.0,
        'k9': 0.1, 'k10': 0.02, 'k11': 0.045, 'k12': 2.27, 'k13': 0.004,
        'gf': 1.0, 'k_epi_on': 0.075, 'Epi_deg': 0.075/0.15, 'k_sil_strength': 1, 
        'K_mybl2_cdc': 0.172, 'h_mybl2_cdc': 4, 'E_sil_thresh': 1.5,
        'k_ecad_prod': 1.2,        
        'k_ecad_deg_m': 0.12,      
        'k_ecad_trans': 0.6,        
        'k_ecad_deg_p': 0.08,        
        'Ecad_ref_scale': 75.0,
        'hills_coeff_silencing': 5,
        'hills_coeff_fzr1_surge': 5,
        'Fzr1_Ecad_thres': 1.5,
        'upreg_k3m_fzr1': 2,
        'Ecad_upreg': 1.73
    }
        
    # --- Create Final Merged Parameter Dictionary ---
    params = {**params_model1_hr, **params_model2_hr}
    return params

def get_initial_conditions():
    y0_yap = [
        0.106, 0.91, 0.429, 0.571, 0.684, 0.316, 0.651, 0.149, 1.044, 3.956, 
        1.391, 0.609, 1.019, 0.981, 310.139, 189.861, 0.013, 0.105, 105.085, 
        3394.915, 6.007, 0.493, 0.622, 2.943, 1.669, 0.196, 0]
    y0_cc = [0.0] * 11 + [1.0, 1.0]
    y0 = y0_yap + y0_cc
    return y0

def take_parameters(params, custom_fak=None):
    # add/substract a 10% variation to each parameter
    varied_params = {}
    for key, value in params.items():
        if key != 'fak_upreg_K':
            variation = 0.01 * value
            varied_value = np.random.normal(value, variation)
            varied_params[key] = varied_value
            continue
            
        if custom_fak is not None:
            varied_params[key] = custom_fak
        else:
            value = params[key]
            variation = 1/12
            varied_value = value + np.random.exponential(variation)
            varied_params[key] = varied_value
    return varied_params

# =============================================================================
# === 3. SIMULATION EXECUTION FUNCTION
# =============================================================================

def run_stiffness_switch_simulation(y0, params, t_start_h, t_switch_h, t_final_h, E_initial, E_final):
    y_initial = y0.copy()
    solutions_all = [] 

    counter_division_events_phase1 = 0
    counter_division_events_phase2 = 0
    
    # --- Phase 1: High Stiffness ---
    E_current = E_initial
    t_end_phase1 = t_switch_h
    t_current_start = t_start_h
    
    while t_current_start < t_end_phase1:
        sol = solve_ivp(
            fun=yap_cell_cycle_model,
            t_span=[t_current_start, t_end_phase1],
            y0=y_initial,
            args=(params, E_current),
            dense_output=True,
            events=division_event,
            method='BDF',
            atol=1e-8, rtol=1e-8
        )
        solutions_all.append(sol)

        if not sol.t_events[0].size:
            t_current_start = t_end_phase1
            y_initial = sol.y[:, -1]
            break 
        
        t_current_start = sol.t_events[0][0]
        y_at_event = sol.y_events[0][0]

        y_initial = y_at_event.copy()
        # Apply division: halve elements according to original logic framework
        for idx in range(len(y_at_event)):
            if idx not in [37, 38, 39] and idx > 26:  # EpiSil index bounds
                y_initial[idx] = y_at_event[idx] / 2.0
        counter_division_events_phase1 += 1

    # --- Phase 2: Low Stiffness ---
    E_current = E_final
    t_end_phase2 = t_final_h

    while t_current_start < t_end_phase2:
        sol = solve_ivp(
            fun=yap_cell_cycle_model,
            t_span=[t_current_start, t_end_phase2],
            y0=y_initial,
            args=(params, E_current),
            dense_output=True,
            events=division_event,
            method='BDF',
            atol=1e-8, rtol=1e-8
        )
        solutions_all.append(sol)

        if not sol.t_events[0].size:
            break 
        
        t_current_start = sol.t_events[0][0]
        y_at_event = sol.y_events[0][0]

        y_initial = y_at_event.copy()
        for idx in range(len(y_at_event)):
            if idx not in [37, 38, 39] and idx > 26:  
                y_initial[idx] = y_at_event[idx] / 2.0
        counter_division_events_phase2 += 1

    # --- Consolidate Results ---
    t_plot_h = np.linspace(t_start_h, t_final_h, 5000) 
    y_plot = np.zeros((len(y0), len(t_plot_h)))
    
    current_t_start_h = t_start_h
    for sol in solutions_all:
        segment_end_time_h = sol.t[-1]
        indices = (t_plot_h >= current_t_start_h) & (t_plot_h <= segment_end_time_h)
        
        if np.any(indices):
            y_plot[:, indices] = sol.sol(t_plot_h[indices])
        
        current_t_start_h = segment_end_time_h
        
    return t_plot_h, y_plot, counter_division_events_phase1, counter_division_events_phase2

# =============================================================================
# === 4. PLOTTING FUNCTION
# =============================================================================

def plot_time_course_for_fak(t_plot_h, y_plot, t_switch_h, p, fak_val):
    """
    Generates time course plots for a specific FAK configuration.
    """
    print(f"Plotting results for fak_upreg_K = {fak_val}...")
    
    fig, axs = plt.subplots(4, 1, figsize=(14, 20), sharex=True)
    
    idx_YAP_nuc = 22
    idx_mybl2 = 24
    idx_yap_mybl2 = 25
    idx_CycB = 28
    idx_C20m = 29
    idx_Cdc20A = 34
    idx_Cdh1 = 33
    idx_EpiSil = 37
    idx_cdh1m = 31
    idx_FAK = 0
    
    # Plot 1: YAP/TAZ and MYBL2
    ratio_YAP = y_plot[idx_YAP_nuc, :] / (y_plot[16, :] + y_plot[17, :])
    axs[0].plot(t_plot_h, ratio_YAP, label='Nuclear/Cytoplasmic YAP Ratio (a.u.)', color='green')
    axs[0].axvline(x=t_switch_h, color='black', linestyle=':', label=f'Stiffness Drop')
    axs[0].set_ylabel('YAP Ratio (a.u.)')
    axs[0].set_title(f'YAP/TAZ Activity (fak_upreg_K = {fak_val})')
    axs[0].legend()
    axs[0].set_ylim(bottom=0)

    # Plot 2: Cell Cycle Oscillators
    axs[1].plot(t_plot_h, y_plot[idx_CycB, :], label='CycB (a.u.)', color='blue')
    axs[1].plot(t_plot_h, y_plot[idx_Cdh1, :], label='Cdh1 (a.u.)', color='purple', linestyle='--')
    axs[1].plot(t_plot_h, y_plot[idx_Cdc20A, :], label='Cdc20A (Active, a.u.)', color='black', linestyle=':') 
    axs[1].axvline(x=t_switch_h, color='black', linestyle=':')
    axs[1].set_ylabel('Protein Conc. (a.u.)')
    axs[1].set_title('Cell Cycle Dynamics (CycB and Cdh1)')
    axs[1].legend()
    axs[1].set_ylim(bottom=0)

    # Plot 3: The Link (Output)
    ax3_twin = axs[2].twinx()
    p1, = axs[2].plot(t_plot_h, y_plot[idx_cdh1m, :], label='Cdh1m (a.u.)', color='crimson')
    p3, = axs[2].plot(t_plot_h, y_plot[idx_cdh1m+1, :], label='Cdh1t (a.u.)', color='purple', linestyle='--')
    axs[2].set_ylabel('Cdh1m / Cdh1t Conc. (a.u.)', color='crimson') 
    axs[2].tick_params(axis='y', labelcolor='crimson')
    axs[2].set_ylim(bottom=0)
    
    p2, = ax3_twin.plot(t_plot_h, y_plot[idx_C20m, :], label='C20m (Cdc20 mRNA, a.u.)', color='orange', linestyle=':')
    ax3_twin.set_ylabel('C20m Conc. (a.u.)', color='orange')
    ax3_twin.tick_params(axis='y', labelcolor='orange')
    ax3_twin.set_ylim(bottom=0)
    
    axs[2].axvline(x=t_switch_h, color='black', linestyle=':')
    axs[2].set_title('The Link: Cdh1 mRNA and C20 mRNA Production') 
    axs[2].legend(handles=[p1, p3, p2]) 

    # Plot 4: Epigenetic Silencing Control
    axs[3].plot(t_plot_h, y_plot[idx_EpiSil, :], label='EpiSil (Cdh1 Silencing, a.u.)', color='black')
    axs[3].plot(t_plot_h, 1/(1+(p['k_sil_strength'] *y_plot[idx_EpiSil, :])), label='Ecad_m (mRNA, a.u.)', color='blue', linestyle='--')
    axs[3].axvline(x=t_switch_h, color='black', linestyle=':')
    axs[3].set_title('Control: EpiSilencing Accumulates Only at Low Stiffness')
    axs[3].set_xlabel('Time (hours)')
    axs[3].set_ylabel('EpiSil Level (a.u.)')
    axs[3].legend()
    axs[3].set_ylim(bottom=0)
    axs[3].set_xlim(left=0) 

    days_ticks = np.arange(0, t_plot_h[-1]/24 + 1, 5)
    days_labels = [f"{int(day)}d" for day in days_ticks]
    plt.xticks(days_ticks * 24, days_labels)
    plt.tight_layout()
    plt.show()
    plt.close()

# =============================================================================
# === 5. MAIN EXECUTION BLOCK (TIME COURSES)
# =============================================================================

if __name__ == '__main__':
    np.random.seed(32)
    
    base_params = get_default_parameters()
    y0 = get_initial_conditions()
    
    E_initial = 308   
    E_final = 1       
    
    t_start_h = 0.0
    t_switch_h = 15 * 24.0   
    t_final_h = 15 * 24.0 + 25 * 24.0 
    
    fak_values_to_test = [1.04, 1.25, 1.4] 
    
    for fak_val in fak_values_to_test:
        print(f"\n--- Running simulation for fak_upreg_K = {fak_val} ---")
        
        current_params = take_parameters(base_params, custom_fak=fak_val)
        
        t_plot, y_plot, div_phase1, div_phase2 = run_stiffness_switch_simulation(
            y0=y0, 
            params=current_params, 
            t_start_h=t_start_h, 
            t_switch_h=t_switch_h, 
            t_final_h=t_final_h, 
            E_initial=E_initial, 
            E_final=E_final
        )
        
        print(f"Division events in phase 1 (Stiff): {div_phase1}")
        print(f"Division events in phase 2 (Soft): {div_phase2}")
        
        plot_time_course_for_fak(t_plot, y_plot, t_switch_h, current_params, fak_val)
# %%
# %%
# =============================================================================
# === SENSITIVITY ANALYSIS ACROSS 3 MODELS
# =============================================================================
def calculate_cell_cycle_metrics(params_in, E_level):
    """
    Runs the simulation for specific parameters on a fixed stiffness
    and returns the stable cell cycle duration (period) in hours
    and the stable nuclear YAP ratio.
    """
    params = params_in.copy()
    y_initial = np.array(get_initial_conditions(), dtype=float)
    
    t_start = 0.0
    t_max_search = 150.0 
    division_times = []
    
    while t_start < t_max_search:
        sol = solve_ivp(
            fun=yap_cell_cycle_model, 
            t_span=[t_start, t_max_search],
            y0=y_initial, 
            args=(params, E_level),
            events=division_event, 
            method='BDF', 
            atol=1e-6, 
            rtol=1e-6
        )
        
        if sol.status == 1: 
            t_event = sol.t_events[0][0]
            division_times.append(t_event)
            
            y_at_event = sol.y_events[0][0]
            y_initial = np.array(y_at_event, dtype=float)
            for idx in range(len(y_at_event)):
                if idx > 26 and idx not in [37, 38, 39]:  
                    y_initial[idx] = y_at_event[idx] / 2.0
            
            t_start = t_event
            
            if len(division_times) >= 5:
                break
        else:
            break

    if len(division_times) >= 2:
        stable_period = division_times[-1] - division_times[-2]
        y_final = y_at_event
    else:
        stable_period = np.nan
        y_final = sol.y[:, -1]
        
    # Calculate YAP ratio (Nuclear / (Phosphorylated + Unphosphorylated Cytoplasmic))
    yap_ratio = y_final[22] / (y_final[16] + y_final[17])
    
    return stable_period, yap_ratio

if __name__ == '__main__':
    np.random.seed(32)
    fak_k_values = np.arange(1.0, 1.5, 0.01)
    
    base_params = get_default_parameters()
    
    models = [
        ("No Fzr1 Surge", 'red', '--', lambda p: {**p, 'model_type': 1}),
        ("No Ecad Recovery", 'blue', '-.', lambda p: {**p, 'model_type': 2}),
        ("Base Model", 'green', '-', lambda p: {**p, 'model_type': 3})
    ]
    
    results_duration = {name: [] for name, _, _, _ in models}
    results_yap = {name: [] for name, _, _, _ in models}

    print("Running parameter sweep over fak_upreg_K for three models...")
    for val in fak_k_values:
        for name, color, style, modifier in models:
            p_set = modifier(base_params)
            p_set['fak_upreg_K'] = val
            
            duration, yap_ratio = calculate_cell_cycle_metrics(p_set, E_level=1)
            results_duration[name].append(duration)
            results_yap[name].append(yap_ratio)
            
        print(f"fak_upreg_K = {val:.2f} processed.")

    # Plotting results in two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    for name, color, style, _ in models:
        # Filter out nan values specifically for the duration plot lines/markers
        valid_indices = ~np.isnan(results_duration[name])
        valid_fak = fak_k_values[valid_indices]
        valid_duration = np.array(results_duration[name])[valid_indices]
        
        # Panel 1: Cell Cycle Duration (Only valid points)
        ax1.scatter(valid_fak, valid_duration, color=color)
        ax1.plot(valid_fak, valid_duration, color=color, linestyle=style, linewidth=2, label=name)
        
        # Panel 2: YAP Ratio (All points)
        ax2.scatter(fak_k_values, results_yap[name], color=color)
        ax2.plot(fak_k_values, results_yap[name], color=color, linestyle=style, linewidth=2, label=name)
    
    # Format top panel
    ax1.set_title('Sensitivity Analysis: Impact of FAK Upregulation Capacity', fontsize=16, pad=15)
    ax1.set_ylabel('Cell Cycle Duration (Hours)', fontsize=14)
    ax1.legend()
    ax1.grid(alpha=0.3, linestyle='--')
    
    # Format bottom panel
    ax2.set_xlabel('FAK Carrying Capacity Parameter ($fak\_upreg\_K$)', fontsize=14)
    ax2.set_ylabel('Nuclear/Cytoplasmic YAP Ratio (a.u.)', fontsize=14)
    ax2.set_xlim([1.0, 1.5])
    ax2.grid(alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.show()
# %%

#=========================================================================================
# Single parameter sensitivity analysis
#=========================================================================================

import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

# =============================================================================
# === 4. SENSITIVITY ANALYSIS 
# =============================================================================

def evaluate_metrics(params, t_start=0, t_switch=15*24, t_final=40*24, E_init=308, E_fin=1):
    """
    Runs the simulation for a given parameter set and returns:
    1. Average cell cycle duration in Phase 2 (after stiffness switch)
    2. Nuclear YAP ratio at the end of the simulation
    """
    y_curr = get_initial_conditions()
    
    # --- Phase 1: High Stiffness ---
    t_curr = t_start
    while t_curr < t_switch:
        sol = solve_ivp(
            fun=yap_cell_cycle_model, t_span=[t_curr, t_switch], y0=y_curr,
            args=(params, E_init), dense_output=False, events=division_event,
            method='BDF', atol=1e-8, rtol=1e-8
        )
        if not sol.t_events[0].size:
            y_curr = sol.y[:, -1]
            break
        
        t_curr = sol.t_events[0][0]
        y_curr = sol.y_events[0][0].copy()
        for idx in range(27, len(y_curr)):
            if idx not in [37, 38, 39]:
                y_curr[idx] /= 2.0

    # --- Phase 2: Low Stiffness ---
    t_curr = t_switch
    div_times = []
    while t_curr < t_final:
        sol = solve_ivp(
            fun=yap_cell_cycle_model, t_span=[t_curr, t_final], y0=y_curr,
            args=(params, E_fin), dense_output=False, events=division_event,
            method='BDF', atol=1e-8, rtol=1e-8
        )
        if not sol.t_events[0].size:
            y_curr = sol.y[:, -1]
            break
        
        t_curr = sol.t_events[0][0]
        div_times.append(t_curr)
        y_curr = sol.y_events[0][0].copy()
        for idx in range(27, len(y_curr)):
            if idx not in [37, 38, 39]:
                y_curr[idx] /= 2.0

    # Calculate average cell cycle duration (Phase 2)
    if len(div_times) > 1:
        avg_cc = np.mean(np.diff(div_times))
    elif len(div_times) == 1:
        avg_cc = div_times[0] - t_switch
    else:
        avg_cc = np.nan

    # Calculate Nuclear YAP Ratio at the end
    # YAPTAZ_p (16), YAPTAZ (17), YAPTAZ_nuc (22), yap_mybl2 (25)
    nuc_yap = y_curr[22] + y_curr[25] 
    total_yap = y_curr[16] + y_curr[17] + y_curr[22] + y_curr[25]
    yap_ratio = nuc_yap / total_yap if total_yap > 0 else np.nan

    return avg_cc, yap_ratio

def run_sensitivity_analysis(variation=0.10):
    base_params = get_default_parameters()
    
    print("Evaluating Base Model...")
    base_cc, base_yap = evaluate_metrics(base_params)
    
    print("\n=== Base Values ===")
    print(f"Base Cell Cycle Duration (Phase 2) : {base_cc:.2f} hours" if not np.isnan(base_cc) else "Base Cell Cycle Duration (Phase 2) : NaN")
    print(f"Base Nuclear YAP Ratio (End of Sim): {base_yap:.4f}\n")

    results = []
    total_params = len([k for k in base_params.keys() if k != 'model_type'])
    print(f"Starting Sensitivity Analysis (+/- {variation*100:.0f}%) for {total_params} parameters...")

    for i, (param_name, base_val) in enumerate(base_params.items()):
        if param_name == 'model_type':
            continue
            
        print(f"Processing ({i+1}/{total_params}): {param_name}")
        
        # +% Change
        params_up = base_params.copy()
        params_up[param_name] = base_val * (1 + variation)
        cc_up, yap_up = evaluate_metrics(params_up)
        
        # -% Change
        params_dn = base_params.copy()
        params_dn[param_name] = base_val * (1 - variation)
        cc_dn, yap_dn = evaluate_metrics(params_dn)
        
        # Calculate % changes
        pct_cc_up = ((cc_up - base_cc) / base_cc * 100) if not np.isnan(cc_up) and not np.isnan(base_cc) else np.nan
        pct_cc_dn = ((cc_dn - base_cc) / base_cc * 100) if not np.isnan(cc_dn) and not np.isnan(base_cc) else np.nan
        
        pct_yap_up = ((yap_up - base_yap) / base_yap * 100) if not np.isnan(yap_up) and not np.isnan(base_yap) else np.nan
        pct_yap_dn = ((yap_dn - base_yap) / base_yap * 100) if not np.isnan(yap_dn) and not np.isnan(base_yap) else np.nan
        
        results.append({
            'Parameter': param_name,
            'CC +10% Change (%)': pct_cc_up,
            'CC -10% Change (%)': pct_cc_dn,
            'YAP +10% Change (%)': pct_yap_up,
            'YAP -10% Change (%)': pct_yap_dn
        })

    df_results = pd.DataFrame(results)

    format_dict = {
        'CC +10% Change (%)': lambda x: f"{x:+.2f}" if pd.notnull(x) else "NaN",
        'CC -10% Change (%)': lambda x: f"{x:+.2f}" if pd.notnull(x) else "NaN",
        'YAP +10% Change (%)': lambda x: f"{x:+.2f}" if pd.notnull(x) else "NaN",
        'YAP -10% Change (%)': lambda x: f"{x:+.2f}" if pd.notnull(x) else "NaN"
    }

    print("\n" + "="*80)
    print(f"{'SENSITIVITY ANALYSIS RESULTS':^80}")
    print("="*80)
    print(df_results.to_string(formatters=format_dict, index=False))
    print("="*80)
    
    return df_results

if __name__ == "__main__":
    df_sensitivity = run_sensitivity_analysis(variation=0.10)
# %%
