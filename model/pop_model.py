#%%
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.size': 16, 'font.family': 'Arial'})

# =============================================================================
# === 1. MODEL DEFINITION (ODE FUNCTION AND EVENTS)
# =============================================================================

def yap_cell_cycle_model(t, y, p, E):
    """
    MERGED MODEL (Time unit: hours)
    
    All parameters are pre-scaled to /h in the parameter dictionary.
    Both Model 1 (YAP/TAZ) and Model 2 (Cell Cycle) derivatives are calculated 
    directly in /h and used as-is.
    """
    
    # --- Unpack All 39 Species ---
    FAK, pFAK, RhoA_GDP, RhoA_GTP, ROCK, ROCK_A, mDia, mDia_A, Myo, Myo_A, \
    LIMK, LIMK_A, Cofilin_p, Cofilin_NP, G_actin, F_actin, YAPTAZ_p, YAPTAZ, \
    LaminA_p, LaminA, NPC, NPC_A, YAPTAZ_nuc, Importin, mybl2, yap_mybl2, \
    SF_remod = y[:27]
    
    Bm, CycB, C20m, Cdc20t, cdhm, cdht, Cdh1, Cdc20A, IEP, cdt1, EpiSil, Ecad_m, Ecad = y[27:]
    
    # =================================================================
    # --- MODULE 1: YAP/TAZ Mechanotransduction ---
    # =================================================================

    basal_activation = p['k_f'] * FAK 
    stiffness_activation = p['k_sf'] * (E / (p['C'] + E)) * FAK 
    R1 = basal_activation + stiffness_activation 
    R2 = p['k_df'] * pFAK 
    
    rho_activation_rate = p['k_fk_rho'] * (p['gamma'] * pFAK**p['n'] + 1) 
    mybl2_feedback_rate = p['k_rho_mybl2'] * yap_mybl2  
    R3 = (rho_activation_rate + mybl2_feedback_rate) * RhoA_GDP 
    R4 = p['k_d_rho'] * RhoA_GTP 

    R5 = p['k_drock'] * ROCK_A 
    R6 = p['k_r_rho'] * RhoA_GTP * ROCK 
    T_ROCKA = (np.tanh(p['sc1'] * (ROCK_A - p['ROCKb'])) + 1) * ROCK_A / 2 

    R7 = p['k_dmdia'] * mDia_A 
    R8 = p['k_m_rho'] * RhoA_GTP * mDia 

    R9 = p['k_mr'] * (p['epsilon'] * T_ROCKA + 1) * Myo - p['k_dmy'] * Myo_A 
    R10 = p['k_lr'] * (p['tau'] * T_ROCKA + 1) * LIMK - p['k_dl'] * LIMK_A 
    R11 = p['k_turn_over'] * Cofilin_p - (p['k_catcofilin'] * LIMK_A * Cofilin_NP) / (p['k_mcofilin'] + Cofilin_NP) 
    R12 = p['k_ra'] * (p['alpha'] * T_ROCKA + 1) * G_actin - (p['k_dep'] + p['k_fc1'] * Cofilin_NP) * F_actin 
    
    R13 = (p['k_CN'] + p['k_CY'] * F_actin * Myo_A) * YAPTAZ_p - p['k_NC'] * YAPTAZ 
    E_cytosol = p['p'] * F_actin**2.6 
    R14 = p['k_fl'] * (E_cytosol / (p['C_Lamin'] + E_cytosol)) * LaminA_p - p['k_rl'] * LaminA 
    R15 = p['k_fNPC'] * LaminA * F_actin * Myo_A * NPC - p['k_rNPCA'] * NPC_A 
    R16 = (p['k_inb_max'] * (Importin/3) + p['k_in'] * NPC_A) * YAPTAZ - p['k_out'] * YAPTAZ_nuc  
    
    production_imp = p['basal_imp_prod'] + p['k_prod_imp'] * (YAPTAZ_nuc**p['h_imp'] / (p['K_imp']**p['h_imp'] + YAPTAZ_nuc**p['h_imp'])) 
    degradation_imp = p['k_deg_imp'] * Importin 
    dImportin_dt_h = production_imp - degradation_imp 
    
    production_mybl2 = p['basal_mybl2_prod'] + p['k_prod_mybl2'] * (YAPTAZ_nuc**p['h_mybl2'] / (p['K_mybl2']**p['h_mybl2'] + YAPTAZ_nuc**p['h_mybl2'])) 
    degradation_mybl2 = p['k_deg_mybl2'] * mybl2 
    net_prod_mybl2 = production_mybl2 - degradation_mybl2 

    d_yap_mybl2_dt_h = p['comp_f_ym'] * mybl2 * YAPTAZ_nuc - p['k_diss_ym'] * yap_mybl2 

    k_remod_f_eff = p['k_remod_f'] if E < p['E_sil_thresh'] else 0.0
    k_remod_d_eff = p['k_remod_d'] if E < p['E_sil_thresh'] else 0.0
    dSF_remod_dt_h = k_remod_f_eff * F_actin - k_remod_d_eff * SF_remod

    production_FAK = p['FAK_prod_rate']

    dFAK_dt_h = R2 - R1 + production_FAK*(1 - (FAK+pFAK)/(p['fak_upreg_K']))* SF_remod - p['deg_FAK']*FAK* SF_remod 
    dpFAK_dt_h = R1 - R2 
    dRhoA_GDP_dt_h = -R3 + R4 
    dRhoA_GTP_dt_h = R3 - R4
    dROCK_dt_h = R5 - R6 
    dROCK_A_dt_h = -R5 + R6
    dmDia_dt_h = R7 - R8 
    dmDia_A_dt_h = -R7 + R8
    dMyo_dt_h = -R9 
    dMyo_A_dt_h = R9
    dLIMK_dt_h = -R10 
    dLIMK_A_dt_h = R10
    dCofilin_p_dt_h = -R11 
    dCofilin_NP_dt_h = R11
    dG_actin_dt_h = -R12 
    dF_actin_dt_h = R12
    dYAPTAZ_p_dt_h = -R13 
    dYAPTAZ_dt_h = R13 - R16
    dLaminA_p_dt_h = -R14 
    dLaminA_dt_h = R14
    dNPC_dt_h = -R15 
    dNPC_A_dt_h = R15
    dYAPTAZ_nuc_dt_h = R16 - d_yap_mybl2_dt_h 
    dmybl2_dt_h = net_prod_mybl2 - d_yap_mybl2_dt_h 
    dyap_mybl2_dt_h = d_yap_mybl2_dt_h 

    derivs_model1_h = [
        dFAK_dt_h, dpFAK_dt_h, dRhoA_GDP_dt_h, dRhoA_GTP_dt_h, dROCK_dt_h, dROCK_A_dt_h,
        dmDia_dt_h, dmDia_A_dt_h, dMyo_dt_h, dMyo_A_dt_h, dLIMK_dt_h, dLIMK_A_dt_h,
        dCofilin_p_dt_h, dCofilin_NP_dt_h, dG_actin_dt_h, dF_actin_dt_h, dYAPTAZ_p_dt_h,
        dYAPTAZ_dt_h, dLaminA_p_dt_h, dLaminA_dt_h, dNPC_dt_h, dNPC_A_dt_h,
        dYAPTAZ_nuc_dt_h, dImportin_dt_h, dmybl2_dt_h, dyap_mybl2_dt_h, dSF_remod_dt_h
    ]

    # =================================================================
    # --- MODULE 2: Cell Cycle ---
    # =================================================================

    cdht_minus_cdh1 = max(0, cdht - Cdh1)
    denom_cdh1 = (p['J3']) + cdht_minus_cdh1
    if denom_cdh1 == 0: denom_cdh1 = 1e-9

    cdc20t_minus_cdc20a = max(0, Cdc20t - Cdc20A)
    denom_cdc20a = (p['J7']) + cdc20t_minus_cdc20a
    if denom_cdc20a == 0: denom_cdc20a = 1e-9

    dBm_dt_h = (((p['k1m'] * p['d']) * p['GF']) / (p['kmm'] + (p['keff'] * p['GF'])) - (p['k1dm'] * p['d']) * Bm) 
    dCycB_dt_h = (p['k1'] * Bm - p['k2a'] * CycB - p['k2b'] * CycB * Cdh1) * p['d'] 

    cycb_hill = (CycB ** p['n_cc']) / ((p['J5']) ** p['n_cc'] + (CycB ** p['n_cc'])) 
    mybl2_activation = (yap_mybl2**p['h_mybl2_cdc']) / (p['K_mybl2_cdc']**p['h_mybl2_cdc'] + yap_mybl2**p['h_mybl2_cdc']) 
    
    basal_prod_c20m = (p['k5am']) 
    cycb_prod_c20m = (p['k5bm'] * cycb_hill) / (p['k5cm'] + (p['GF'] * p['j5c']))
    total_prod_c20m = (basal_prod_c20m + cycb_prod_c20m) * mybl2_activation 
    dC20m_dt_h = (total_prod_c20m - p['k5dm'] * C20m) * p['d'] 

    dCdc20t_dt_h = (p['k5a'] * C20m - p['k6'] * Cdc20t) * p['d'] 

    silencing_factor = 1.0 / (1.0 + (SF_remod / p['k_sil_strength'])**p['hills_coeff_silencing'])
    ecad_normalized = max(0.0, (Ecad / p['Ecad_ref_scale']) - 1.0)
    fzr1_hill_surge = (ecad_normalized**p['hills_coeff_fzr1_surge']) / (ecad_normalized**p['hills_coeff_fzr1_surge'] + (p['Fzr1_Ecad_thres'])**p['hills_coeff_fzr1_surge'])

    model_type = p.get('model_type', 3)
    if model_type == 1:
        fzr1_hill_surge = 0.0
    elif model_type == 2:
        silencing_factor = 1.0

    k3m_effective = p['k3m'] * (1.0 + p['upreg_k3m_fzr1'] * fzr1_hill_surge) if E < p['E_sil_thresh'] else p['k3m']

    dcdhm_dt_h = (k3m_effective - p['k3dm'] * cdhm) * p['d'] 
    dcdht_dt_h = (p['k3a'] * cdhm - p['k3dt'] * cdht) * p['d'] 

    ke_effective = p['k3'] 
    term1_cdh1 = ((ke_effective) + (p['k3b'] * Cdc20A)) * cdht_minus_cdh1 / denom_cdh1 
    term2_cdh1 = p['k4'] * CycB * Cdh1 / ((p['J4']) + Cdh1) 
    dCdh1_dt_h = (term1_cdh1 - term2_cdh1 - p['k3dt'] * Cdh1) * p['d'] 

    term1_cdc20a = p['k7'] * IEP * cdc20t_minus_cdc20a / denom_cdc20a
    term2_cdc20a = p['k8'] * p['Mad'] * Cdc20A / ((p['J8']) + Cdc20A)
    dCdc20A_dt_h = (term1_cdc20a - term2_cdc20a - p['k6'] * Cdc20A) * p['d'] 

    dIEP_dt_h = (p['k9'] * (CycB) * (1 - IEP) - p['k10'] * IEP) * p['d'] 
    dcdt1_dt_h = (p['k11'] - p['k12'] * CycB * cdt1 - p['k13'] * cdt1) * p['d'] 

    dEpiSil_dt_h = p['k_epi_on'] * ((SF_remod**3)/(p['threshold_remod_to_DNMT_act']**3 + SF_remod**3)) - p['Epi_deg'] * EpiSil #this equation and associated parameters have not been used in the model

    k_surge_ecad_prod = (p['k_ecad_prod'] + ((p['Ecad_upreg']*p['k_ecad_prod'])*silencing_factor)) if E < p['E_sil_thresh'] else p['k_ecad_prod']
    dEcad_m_dt_h = (k_surge_ecad_prod  - p['k_ecad_deg_m'] * Ecad_m) * p['d']
    dEcad_dt_h = (p['k_ecad_trans'] * Ecad_m - p['k_ecad_deg_p'] * Ecad) * p['d']

    derivs_model2_h = [
        dBm_dt_h, dCycB_dt_h, dC20m_dt_h, dCdc20t_dt_h, dcdhm_dt_h, dcdht_dt_h,
        dCdh1_dt_h, dCdc20A_dt_h, dIEP_dt_h, dcdt1_dt_h, dEpiSil_dt_h, dEcad_m_dt_h, dEcad_dt_h
    ]
    
    return derivs_model1_h + derivs_model2_h

def division_event(t, y, p, E):
    return y[28] - 0.1 

division_event.terminal = True 
division_event.direction = -1 

# =============================================================================
# === 2. SETUP FUNCTIONS
# =============================================================================

def get_default_parameters():
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
        'k_remod_f': 1e-7*3600, 'k_remod_d': 1e-5*3600, 'FAK_prod_rate': 0.35*3600, 'deg_FAK': 0.15*3600, 'fak_upreg_K': 1, 'threshold_remod_to_DNMT_act': 0.5
    }
    
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

def take_parameters(params):
    varied_params = {}
    for key, value in params.items():
        if key == 'model_type':
            varied_params[key] = value
            continue
        if key != 'fak_upreg_K':
            variation = 0.0001 * value
            varied_value = np.random.normal(value, variation)
            varied_params[key] = varied_value
            continue
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
        for idx in range(len(y_at_event)):
            if idx not in [37, 38, 39] and idx > 26:  
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
#%%
# =============================================================================
# === 4. EXECUTE PIPELINE FOR MODELS 1, 2, AND 3
# =============================================================================

y0 = get_initial_conditions()
E_initial = 308 
E_final = 1     
t_start_h = 0.0
t_switch_h = 5 * 24.0   
t_final_h = 5 * 24.0 + 10 * 24.0 

# Dictionary to store weighted average birth rates for combined plotting
all_avg_division_rates = {}

for model_type in [3,2,1]:
    np.random.seed(42)  # Reset the random seed for reproducibility
    print(f"\n============================================================")
    print(f"RUNNING ANALYSIS FOR MODEL TYPE {model_type}")
    print(f"============================================================")

    array_fak_upreg_K = []
    fak_levels_array = []
    nuclear_yap_array = []
    cdh1_levels_array = []
    avg_divisions_phase1 = []
    avg_divisions_phase2 = []

    for _i in range(150):
        print(f"--- Model {model_type} | Sample {_i + 1}/250 ---")
        params = get_default_parameters()
        params['model_type'] = model_type
        varied_params = take_parameters(params)
        
        fak_upreg_K_values = varied_params['fak_upreg_K']
        array_fak_upreg_K.append(fak_upreg_K_values)

        t_plot, y_plot, counter_division_events_phase1, counter_division_events_phase2 = run_stiffness_switch_simulation(
            y0=y0, 
            params=varied_params,
            t_start_h=t_start_h,
            t_switch_h=t_switch_h,
            t_final_h=t_final_h,
            E_initial=E_initial,
            E_final=E_final
        )

        fak_levels = y_plot[0, -1] + y_plot[1, -1]
        fak_levels_array.append(fak_levels)
        ratio_YAP = y_plot[22, :] / (y_plot[16, :] + y_plot[17, :])
        nuclear_yap_level = ratio_YAP[-1]
        nuclear_yap_array.append(nuclear_yap_level)
        cdh_level = y_plot[31, :][-1]
        cdh1_levels_array.append(cdh_level)
        avg_divisions_phase1.append(counter_division_events_phase1)
        avg_divisions_phase2.append(counter_division_events_phase2)

    # --- Plotting Initial Distributions & Scatters ---
    plt.figure()
    sns.histplot(x=array_fak_upreg_K, bins=50, color='skyblue', edgecolor='black')
    plt.xlabel('fak_upreg_K Values', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(f'Distribution of fak_upreg_K (Model {model_type})', fontsize=16)
    plt.show()

    plt.figure()
    sns.scatterplot(y=nuclear_yap_array, x=fak_levels_array, color='black', linewidth=0)
    plt.xlabel('FAK Levels on Soft Substrate', fontsize=14)
    plt.ylabel('Nuclear to Cytoplasmic YAP Ratio', fontsize=14)
    plt.title(f'Model {model_type}', fontsize=16)
    plt.show()

    plt.figure()
    sns.scatterplot(y=cdh1_levels_array, x=fak_levels_array, color='black', linewidth=0)
    plt.xlabel('FAK Levels on Soft Substrate', fontsize=14)
    plt.ylabel('Cdh1 Levels on Soft Substrate', fontsize=14)
    plt.title(f'Model {model_type}', fontsize=16)
    plt.show()

    plt.figure()
    sns.scatterplot(y=[x/10 for x in avg_divisions_phase2], x=nuclear_yap_array, color='black', linewidth=0)
    plt.ylabel('Growth rate (divisions per day)', fontsize=14)
    plt.xlabel('Nuclear to Cytoplasmic YAP Ratio', fontsize=14)
    plt.title(f'Model {model_type}', fontsize=16)
    plt.show()

    sns.set_style("whitegrid", {'axes.grid' : False})
    plt.figure(figsize=(5, 4))
    sns.histplot(x=nuclear_yap_array, bins=50, color='grey', edgecolor='black')
    plt.xlabel('Nuclear to Cytoplasmic YAP Ratio', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title(f'Distribution of nuclear YAP Levels (Model {model_type})', fontsize=16)
    plt.show()

    # --- Clonal Selection Dynamics ---
    # Choose ALL clones, not just the >0 ones
    #arr_for_selection = avg_divisions_phase2.copy()
    arr_for_selection = []
    array_fak_upreg_K_selection = []
    for i,j in enumerate(avg_divisions_phase2):
        if j>0:
            arr_for_selection.append(j)
            array_fak_upreg_K_selection.append(array_fak_upreg_K[i])

    #np.random.shuffle(arr_for_selection)

    if not arr_for_selection:
        print(f"No cell divisions observed for Model {model_type}. Skipping clonal selection plots.\n")
        continue

    DIVISION_RATES = np.array([x/10 for x in arr_for_selection])  
    
    # Sigmoid death rate starting at 5% up to 20%, thresholded at 1.325
    DEATH_FRACTION = 0.20#0.05 + (0.15 / (1.0 + np.exp(-10 * (np.array(array_fak_upreg_K_selection) - 1.325))))

    NUM_CLONES = len(DIVISION_RATES)
    INITIAL_POPULATION_PER_CLONE = 5
    TOTAL_INITIAL_POPULATION = INITIAL_POPULATION_PER_CLONE * NUM_CLONES
    NUM_TIME_STEPS = 90  
    TIME_PERIOD = 1 

    population_history = pd.DataFrame(
        0.0,
        index=range(NUM_TIME_STEPS + 1),
        columns=[f'Clone {i+1} (Rate: {rate:.3f})' for i, rate in enumerate(DIVISION_RATES)]
    )

    initial_population = np.full(NUM_CLONES, INITIAL_POPULATION_PER_CLONE, dtype=float)
    population_history.iloc[0] = initial_population

    for t in range(NUM_TIME_STEPS):
        current_population = population_history.iloc[t].values
        next_population = np.zeros(NUM_CLONES)

        NET_GROWTH_RATES = DIVISION_RATES * (1 - DEATH_FRACTION)
        growth_factors = 1 + NET_GROWTH_RATES
        total_grown_population = current_population * growth_factors

        total_current_growth = np.sum(total_grown_population)
        scaling_factor = TOTAL_INITIAL_POPULATION / total_current_growth if total_current_growth > 0 else 0

        next_population = total_grown_population * scaling_factor
        population_history.iloc[t+1] = next_population

    # Plot 6: Population History
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))

    population_history.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        colormap='rainbow',
        legend=False # Legend removed
    )

    ax.set_title(f'Clonal Selection Dynamics (Model {model_type})', fontsize=16, fontweight='bold')
    ax.set_xlabel(f'Time Step ({TIME_PERIOD}-Day Period)', fontsize=12)
    ax.set_ylabel('Total Cell Population (Constant Carrying Capacity)', fontsize=12)
    ax.set_xticks(np.arange(0, NUM_TIME_STEPS + 1, 5)) 
    ax.set_xticklabels(np.arange(0, (NUM_TIME_STEPS + 1) * TIME_PERIOD, TIME_PERIOD * 5))
    ax.tick_params(axis='x', rotation=0)

    final_composition = population_history.iloc[-1].idxmax()
    textstr = (f'Outcome (Model {model_type}):\n{final_composition} dominates.\n'
               f'Net Growth Rate = Birth Rate * (1 - Variable Death Fraction)\n'
               f'Accounts for {population_history.iloc[-1].max()/TOTAL_INITIAL_POPULATION:.1%} of final pop.')

    ax.text(
        0.2, -0.5, textstr, 
        transform=ax.transAxes, 
        fontsize=12, 
        verticalalignment='top', 
        bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.6)
    )

    plt.tight_layout() 
    plt.show()

    # Track weighted average birth rate over time for the combined plot
    avg_division_rate_over_time = population_history.apply(
        lambda row: np.sum(row.values * NET_GROWTH_RATES) / np.sum(row.values) if np.sum(row.values) > 0 else 0.0, axis=1
    )
    all_avg_division_rates[model_type] = avg_division_rate_over_time

# Plot 7: Combined Weighted Average Birth Rate Over Time
sns.set_style("whitegrid", {'axes.grid' : False})
plt.figure(figsize=(10, 6))

for m_type, avg_rates in all_avg_division_rates.items():
    plt.plot(avg_rates, marker='o', label=f'Model {m_type}')

plt.title('Weighted Average Birth Rate Over Time (All Models)', fontsize=16)
plt.xlabel(f'Time Step ({TIME_PERIOD}-Day Period)', fontsize=12)
plt.ylabel('Average Division Rate (Divisions/10 Days)', fontsize=12)
plt.grid(True)
plt.xticks(np.arange(0, NUM_TIME_STEPS + 1, 5), np.arange(0, (NUM_TIME_STEPS + 1) * TIME_PERIOD, TIME_PERIOD * 5))
plt.legend()
plt.show()
# %%
