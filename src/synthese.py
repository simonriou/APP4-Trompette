import os
import json
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.signal
from scipy.io import wavfile
import matplotlib.pyplot as plt

def charger_parametres(fichier="config.json"):
    with open(fichier, 'r') as f:
        return json.load(f)

def p0_func(x, params_p0):
    return params_p0["amplitude"] * np.exp(-0.5 * ((x - params_p0["x0"]) / params_p0["sigma"])**2)

def dp0_dt_func(x, params_dp0):
    return params_dp0["amplitude"]

def calcul_Zc_kc_eff(omega, r1, r2, params_phys, params_dissip):
    rho0 = params_phys["rho0"]
    c = params_phys["c"]
    gamma = params_phys["gamma"]
    Cp = params_phys["Cp"]
    
    rm = (r1 + r2) / 2.0
    if not params_dissip["activer"]:
        return rho0 * c, omega / c

    if params_dissip["mode_entree"] == "coefficients":
        mu = params_dissip["mu"]
        kappa = params_dissip["kappa"]
    else:
        omega_ref = 2 * np.pi * params_dissip["f_ref"]
        delta_v = params_dissip["delta_v"]
        delta_t = params_dissip["delta_t"]
        mu = (delta_v**2) * rho0 * omega_ref
        kappa = (delta_t**2) * rho0 * Cp * omega_ref

    lv = mu / (rho0 * c)
    Pr = (mu * Cp) / kappa
    
    rv = rm * np.sqrt(omega / (c * lv))
    rt = rv * np.sqrt(Pr)
    sqrt_mj = np.sqrt(-1j)
    
    rho_s = rho0 * (1.0 + 2.0 * sqrt_mj / rv)
    chi_s = 1.0 / (rho0 * c**2)
    chi_t = chi_s * (1.0 + 2.0 * (gamma - 1.0) * sqrt_mj / rt)
    
    k = omega * np.sqrt(rho_s * chi_t)
    rho_c_eff = np.sqrt(rho_s / chi_t)
    return rho_c_eff, k

def matrice_transfert_segment(omega, L, r1, r2, phys, dissip):
    rho_c_eff, k = calcul_Zc_kc_eff(omega, r1, r2, phys, dissip)
    
    if np.abs(r1 - r2) < 1e-9:
        Zc = rho_c_eff / (np.pi * r1**2)
        return np.array([
            [np.cos(k * L), 1j * Zc * np.sin(k * L)],
            [1j * (1.0 / Zc) * np.sin(k * L), np.cos(k * L)]
        ], dtype=complex)
    
    x1 = r1 * L / (r2 - r1)
    x2 = r2 * L / (r2 - r1)
    
    M11 = (r2 / r1) * np.cos(k * L) - np.sin(k * L) / (k * x1)
    M12 = 1j * rho_c_eff / (np.pi * r1 * r2) * np.sin(k * L)
    M21 = (np.pi * r1 * r2 / rho_c_eff) * (
        1j * np.sin(k * L) * (1.0 + 1.0 / (k**2 * x1 * x2)) 
        + (np.cos(k * L) / (1j * k)) * (1.0 / x1 - 1.0 / x2)
    )
    M22 = (r1 / r2) * np.cos(k * L) + np.sin(k * L) / (k * x2)
    
    return np.array([[M11, M12], [M21, M22]], dtype=complex)

def calcul_Zr(omega, r_out, phys, ray):
    if ray["type"] == "aucun":
        return 0j
        
    rho0 = phys["rho0"]
    c = phys["c"]
    S = np.pi * r_out**2
    k0 = omega / c
    negliger_Rr = ray["negliger_Rr"]

    if ray["type"] == "baffle":
        Rr = ((k0 * r_out)**2) / 2.0 if not negliger_Rr else 0.0
        Xr = (8.0 * k0 * r_out) / (3.0 * np.pi)
    else:
        Rr = ((k0 * r_out)**2) / 4.0 if not negliger_Rr else 0.0
        Xr = 0.6133 * k0 * r_out
        
    return (rho0 * c / S) * (Rr + 1j * Xr)

def simulation_et_audio():
    params = charger_parametres("config.json")
    phys = params["physique"]
    c = phys["c"]
    
    z_list = np.array(params["geometrie"]["z_list"])
    r_list = np.array(params["geometrie"]["r_list"])
    L_tot = z_list[-1] - z_list[0]
    r_out = r_list[-1]
    
    n_modes = params["simulation"]["n_modes"]
    fs = params["simulation"]["fs"]
    duree = params["simulation"]["duree"]
    ray = params["rayonnement"]
    dissip = params["dissipation"]
    ecoute = params["ecoute"]
    type_ecoute = ecoute["type"]
    
    dossier_sortie = "output"
    os.makedirs(dossier_sortie, exist_ok=True)
    nom_base = f"perce_cascade_{type_ecoute}_{len(z_list)-1}segments"
    fichier_audio = os.path.join(dossier_sortie, f"{nom_base}.wav")
    fichier_freqs = os.path.join(dossier_sortie, f"{nom_base}_resultats.txt")

    t = np.arange(0, duree, 1.0 / fs)
    p_audio = np.zeros_like(t)
    resultats_modaux = []
    
    print(f"Balayage fréquentiel pour trouver les pics sur {len(z_list)-1} segments...")
    
    f_scan = np.arange(20.0, 2500.0, 0.5) 
    mag_Zin_scan = np.zeros_like(f_scan)
    
    for idx, f in enumerate(f_scan): # calculer la matrice de transfert totale & l'impédance d'entrée
        omega = 2 * np.pi * f
        M_tot = np.eye(2, dtype=complex)
        for i in range(len(z_list) - 1):
            M_tot = M_tot @ matrice_transfert_segment(omega, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
        Zr = calcul_Zr(omega, r_out, phys, ray)
        Zin_val = (M_tot[0, 0] * Zr + M_tot[0, 1]) / (M_tot[1, 0] * Zr + M_tot[1, 1])
        mag_Zin_scan[idx] = np.abs(Zin_val)
        
    # Sécurité anti-NaN au cas où les très basses fréquences causeraient des instabilités
    mag_Zin_scan = np.nan_to_num(mag_Zin_scan, nan=0.0)

    # On cherche efficacement les pics de l'impédance d'entrée pour identifier les fréquences de résonance
    pics_idx, _ = scipy.signal.find_peaks(mag_Zin_scan)
    pics_idx = pics_idx[:n_modes] 
    
    if len(pics_idx) == 0:
        print("Erreur : Aucun pic de résonance trouvé. Le vecteur d'impédance est probablement plat ou rempli de zéros.")
        return

    print(f"{len(pics_idx)} modes identifiés. Calcul des racines complexes...")
    
    for n, idx_pic in enumerate(pics_idx, 1): # On fait la résolution modale en partant des fréquences de résonance identifiées
        f_guess = f_scan[idx_pic]
        omega_0 = 2 * np.pi * f_guess
        
        def systeme_ferme(omega):
            M_tot = np.eye(2, dtype=complex)
            for i in range(len(z_list) - 1):
                M_tot = M_tot @ matrice_transfert_segment(omega, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
            Zr = calcul_Zr(omega, r_out, phys, ray)
            return M_tot[1, 0] * Zr + M_tot[1, 1]
            
        omega_complexe = optimize.newton(systeme_ferme, omega_0 + 0j, tol=1e-6, maxiter=200, disp=False)
        omega_rn = np.real(omega_complexe)
        sigma_n = np.abs(np.imag(omega_complexe))
        
        M_tot_reel = np.eye(2, dtype=complex)
        for i in range(len(z_list) - 1):
            M_tot_reel = M_tot_reel @ matrice_transfert_segment(omega_rn, z_list[i+1]-z_list[i], r_list[i], r_list[i+1], phys, dissip)
        Zr_reel = calcul_Zr(omega_rn, r_out, phys, ray)
        Zin_val = (M_tot_reel[0, 0] * Zr_reel + M_tot_reel[0, 1]) / (M_tot_reel[1, 0] * Zr_reel + M_tot_reel[1, 1])
        mag_Zin = np.abs(Zin_val)
        
        resultats_modaux.append((n, omega_rn / (2 * np.pi), sigma_n, mag_Zin))

        # On calcule les déformées modales en partant de la condition de sortie et en remontant vers l'entrée
        X_frontieres = [np.array([calcul_Zr(omega_complexe, r_out, phys, ray), 1.0], dtype=complex)] # les vecteurs P, U aux frontières (seulement la sortie ici)
        
        for i in range(len(z_list) - 1, 0, -1):
            L_i = z_list[i] - z_list[i-1]
            M_i = matrice_transfert_segment(omega_complexe, L_i, r_list[i-1], r_list[i], phys, dissip)
            X_frontieres.insert(0, M_i @ X_frontieres[0])

        max_P = max(np.abs(X[0]) for X in X_frontieres)
        if max_P > 0:
            X_frontieres = [X / max_P for X in X_frontieres] # normalisation des déformées modales pour éviter les débordements numériques
            
        def r_func(z): # fonction de rayon interpolée pour les segments
            i = min(np.searchsorted(z_list, z, side='right') - 1, len(z_list) - 2)
            i = max(0, i)
            return r_list[i] + (z - z_list[i]) / (z_list[i+1] - z_list[i]) * (r_list[i+1] - r_list[i])
            
        def eval_phi_n_complex(z): # fonction de déformée modale complexe
            i = min(np.searchsorted(z_list, z, side='right') - 1, len(z_list) - 2)
            i = max(0, i)
            if np.isclose(z, z_list[i+1]): return X_frontieres[i+1][0]
            M_part = matrice_transfert_segment(omega_complexe, z_list[i+1] - z, r_func(z), r_list[i+1], phys, dissip)
            return (M_part @ X_frontieres[i+1])[0]

        # Plutôt que de prendre directement la valeur de phi_n à z=0, on corrige la phase pour que le maximum de la déformée modale soit réel et positif
        z_grid = np.linspace(0, L_tot, 50)
        phi_grid = [eval_phi_n_complex(z_val) for z_val in z_grid]
        phase_correction = np.angle(phi_grid[np.argmax(np.abs(phi_grid))])

        def eval_phi_n(z):
            return np.real(eval_phi_n_complex(z) * np.exp(-1j * phase_correction))
            
        N_n, prod_scal_a, prod_scal_b = 0.0, 0.0, 0.0
        
        for i in range(len(z_list) - 1):
            z_a, z_b = z_list[i], z_list[i+1]
            def int_N(z_array): return np.array([np.pi * r_func(z)**2 * eval_phi_n(z)**2 for z in z_array])
            def int_A(z_array): return np.array([np.pi * r_func(z)**2 * p0_func(z, params["excitation"]["p0"]) * eval_phi_n(z) for z in z_array])
            def int_B(z_array): return np.array([np.pi * r_func(z)**2 * dp0_dt_func(z, params["excitation"]["dp0_dt"]) * eval_phi_n(z) for z in z_array])
            
            N_n += integrate.fixed_quad(int_N, z_a, z_b, n=7)[0]
            prod_scal_a += integrate.fixed_quad(int_A, z_a, z_b, n=7)[0]
            prod_scal_b += integrate.fixed_quad(int_B, z_a, z_b, n=7)[0]
            
        if N_n == 0: continue
        a_n = prod_scal_a / N_n
        b_n = prod_scal_b / (omega_rn * N_n)
        
        if type_ecoute == "interne":
            x_obs = ecoute.get("x_obs", 0.3)
            q_n_t = (a_n * np.cos(omega_rn * t) + b_n * np.sin(omega_rn * t)) * np.exp(-sigma_n * t)
            p_audio += eval_phi_n(x_obs) * q_n_t
            
        elif type_ecoute == "rayonne":
            D = 2.0 if ray["type"] == "baffle" else 4.0
            r_obs = ecoute.get("r", 1.0)
            facteur_rad = (np.pi * r_out**2) / (D * np.pi * r_obs)
            
            dz = 1e-5
            dphi_n_dx_L = (eval_phi_n(z_list[-1]) - eval_phi_n(z_list[-1] - dz)) / dz
            
            t_retard = t - r_obs/c
            masque = t_retard >= 0
            t_val = t_retard[masque]
            
            q_n_t_retard = np.zeros_like(t)
            q_n_t_retard[masque] = (a_n * np.cos(omega_rn * t_val) + b_n * np.sin(omega_rn * t_val)) * np.exp(-sigma_n * t_val)
            p_audio += -facteur_rad * dphi_n_dx_L * q_n_t_retard

    max_val = np.max(np.abs(p_audio))
    p_audio_norm = np.int16((p_audio / max_val) * 32767) if max_val > 0 else np.int16(p_audio)
    wavfile.write(fichier_audio, fs, p_audio_norm)
    
    with open(fichier_freqs, 'w') as f:
        f.write("Mode, Frequence(Hz), Taux_Amort(s^-1), Module_Zin(Pa.s.m^-3)\n")
        for n, f_n, sigma_n, mag_Zin in resultats_modaux:
            f.write(f"{n}, {f_n:.3f}, {sigma_n:.3f}, {mag_Zin:.3f}\n")

    print(f"-> Audio généré : {fichier_audio}")
    print(f"-> Résultats sauvegardés : {fichier_freqs}")

    plt.figure(figsize=(10, 4))
    plt.plot(t, p_audio, color='C4')
    plt.title(f"Pression {type_ecoute} - Perce Trompette")
    plt.xlabel("Temps (s)")
    plt.ylabel("Pression (Pa)")
    plt.xlim(0, 0.5)
    plt.grid(True, alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    simulation_et_audio()