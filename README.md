# Pulsed Off-Resonance MT with Intrinsic Blood Suppression - EPG-MT code

This repository contains MATLAB simulation code for the paper:

**"Pulsed Off-Resonance MT with Intrinsic Blood Suppression: Modeling with an Extended EPG Framework and Experimental Results"**  
Authors: Qianxue Shan, Vincent Wong, Ziqiang Yu, Zijian Gao, Qiuyi Shen, Chun Liu, Queenie Chan, Winnie Chu, Weitian Chen

---

## Description

This repository provides the simulation framework for magnetization transfer (MT) imaging with intrinsic blood suppression using a pulsed off-resonance MT approach. The approach is modeled using an Extended Phase Graph (EPG) framework, enabling efficient exploration of sequence and physiological parameters.  

### File List
- `Figure2.m`: Generates results for Figure 2 in the paper.  
  **Description:** Validates the off-resonance EPG-MT simulation against the 1000-isochromat Bloch-McConnell (BM) model using liver and blood parameters. The EPG-MT model produces consistent results with significantly improved computational efficiency.  

- `Figure3.m`: Generates results for Figure 3 in the paper.  
  **Description:** Analyzes blood suppression based on flow direction, gradient amplitude ($ G_z $), pulse number ($ N_p $), and $B_1$ inhomogeneity. Blood suppression depends on the flow direction and increases with greater $G_z$ and $N_p$. $B_1$ inhomogeneity has partial effects.  
  

- `Figure4.m`: Generates results for Figure 4 in the paper.  
  **Description:** Investigates blood suppression across $\omega_1$, $\Delta\omega$, $t_p$, and $t_d$.  
 
### Supporting Functions
- `epgmt_Ex_Relax.m`: Models exchange and relaxation effects in the EPG framework.  
- `epgmt_Flow.m`: Simulates flow-induced phase shifts.  
- `epgmt_Grad.m`: Applies gradient-induced dephasing to the EPG states.  
- `epgmt_RF.m`: Simulates RF pulse effects, including off-resonance and phase cycling.  
- `RF_MT.m`: Computes saturation rates for different lineshapes, including SuperLorentzian, Gaussian, and Lorentzian.

---

## Citation

If you use this code or results in your research, please cite our paper:

> Shan Q, Wong V, Yu Z, Gao Z, Shen Q, Liu C, Chan Q, Chu WCW, Chen W. Pulsed Off-Resonance MT with Intrinsic Blood Suppression: Modeling with an Extended EPG Framework and Experimental Results. ISMRM2026

---

## License

This code is for academic and research purposes only. Commercial use or redistribution is strictly prohibited without explicit written permission from the authors.
