# MPF Mapping with Intrinsic Flow Suppression using Pulsed Saturation

This repository contains MATLAB simulation code for the off-resonance Extended Phase Graph with Magnetization Transfer (EPG-MT) framework, supporting both conference and journal publications.

---

## Publications

### Journal Paper (MRM)
**"Macromolecular Proton Fraction Mapping with Intrinsic Flow Suppression using Pulsed Saturation"**  
*Magnetic Resonance in Medicine*, 2026

This paper presents MPF-SPS (MPF mapping using Short Pulsed Saturation), a method that integrates pulsed off-resonance MT preparation with intrinsic blood signal suppression for quantitative macromolecular proton fraction mapping.

### Conference Abstract (ISMRM 2026)
**"Pulsed Off-Resonance MT with Intrinsic Blood Suppression: Modeling with an Extended EPG Framework and Experimental Results"**  
*ISMRM 2026*

---

## Repository Structure

├── paper/ # Code for MRM journal paper
│ ├── Figure2.m # Validation of MT rotation matrix formalism
│ ├── Figure3.m # Validation of off-resonance EPG-MT model
│ ├── Figure4.m # Flow suppression vs direction, gradient, Np, B1
│ ├── Figure5.m # Flow suppression vs ω1, Δω, tp, td
│ ├── Figure6.m # Tissue-dependent flow suppression
│ ├── Figure7.m # Partial volume effect reduction
│ ├── Figure8.m # Rmpfsps sensitivity to tissue parameters
│ └── ...
├── abstract/ # Code for ISMRM abstract
│ ├── Figure2.m
│ ├── Figure3.m
│ ├── Figure4.m
│ └── ...
├── functions/ # Shared supporting functions
│ ├── epgmt_Ex_Relax.m # Exchange and relaxation in EPG framework
│ ├── epgmt_Flow.m # Flow-induced phase shifts
│ ├── epgmt_Grad.m # Gradient-induced dephasing
│ ├── epgmt_RF.m # RF pulse effects with off-resonance
│ └── RF_MT.m # MT saturation rates for different lineshapes
└── README.md

---

## Paper Simulation Studies

The `paper/` folder contains code reproducing all simulation studies from the MRM paper:

| Figure | Description |
|--------|-------------|
| Figure 2 | Validates the proposed MT rotation matrix formalism (Eq. 1) against full Bloch-McConnell equations across different frequency offsets and RF amplitudes |
| Figure 3 | Validates the off-resonance EPG-MT model and demonstrates flow suppression for white matter and blood at various velocities |
| Figure 4 | Analyzes flow suppression performance as a function of flow direction (θ), gradient strength (Gz), number of pulses (Np), and B1 scale |
| Figure 5 | Investigates flow suppression across ω1, Δω, tp, and td parameter space for sequence optimization |
| Figure 6 | Demonstrates tissue-dependent flow suppression characteristics showing blood attenuates faster than brain parenchyma |
| Figure 7 | Validates reduction of partial volume effects from intravascular signal contamination |
| Figure 8 | Characterizes Rmpfsps sensitivity to tissue parameters (T1a, T2a, T2b, kba, fb) |

---

## Abstract Simulation Studies

The `abstract/` folder contains code from the original ISMRM 2026 submission:

| Figure | Description |
|--------|-------------|
| Figure 2 | Validates EPG-MT simulation against 1000-isochromat Bloch-McConnell model |
| Figure 3 | Analyzes blood suppression based on flow direction, gradient amplitude, pulse number, and B1 inhomogeneity |
| Figure 4 | Investigates blood suppression across ω1, Δω, tp, and td |

---

## Key Features

- **Off-resonance EPG-MT Framework**: Extends the EPG formalism to incorporate off-resonance MT effects using rotation theory, providing efficient simulation compared to isochromat-based methods
- **Intrinsic Flow Suppression**: Achieves blood signal attenuation through cumulative phase dispersion from gradient pulses combined with RF pulse trains
- **Quantitative MPF Mapping**: Enables macromolecular proton fraction quantification insensitive to free water pool parameters (T1a, T2a)

---

## Requirements

- MATLAB R2024a or later (MathWorks, Natick, MA, USA)

---

## Citation

If you use this code in your research, please cite:

**Journal Paper:**
> Author One, Author Two, Author Three. Macromolecular Proton Fraction Mapping with Intrinsic Flow Suppression using Pulsed Saturation. *Magnetic Resonance in Medicine*. 2026.

**Conference Abstract:**
> Shan Q, Wong V, Yu Z, Gao Z, Shen Q, Liu C, Chan Q, Chu WCW, Chen W. Pulsed Off-Resonance MT with Intrinsic Blood Suppression: Modeling with an Extended EPG Framework and Experimental Results. *ISMRM 2026*.

---

## License

This code is for academic and research purposes only. Commercial use or redistribution is strictly prohibited without explicit written permission from the authors.

---

## Contact

For questions regarding this code, please contact:
- Weitian Chen: wtchen@cuhk.edu.hk
- Department of Imaging and Interventional Radiology, The Chinese University of Hong Kong
