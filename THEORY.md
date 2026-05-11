# Theory — Amine-Wash CO₂ Capture

A reader-friendly walkthrough of the physics and chemistry implemented in the digital twin. Read this if you want to understand *why* the model is built the way it is before diving into the code.

---

## 1. Why amine wash?

Of all post-combustion CO₂-capture technologies, **aqueous amine scrubbing** is the only one with multi-decade industrial track record and proven scalability. Originally developed for natural-gas sweetening in the 1930s (Bottoms patent, 1930), it was deployed at megawatt scale for CCS at SaskPower's Boundary Dam (Cansolv DC-103, 2014) and NRG's Petra Nova (KS-1, 2017–2020).

The defining strength of the technology is its **chemical-equilibrium reversibility**: CO₂ reacts with the amine at low temperature (40–60 °C, absorber) and is released at higher temperature (110–125 °C, stripper). This temperature swing is the energy lever that makes regeneration possible — but it is also the dominant cost driver, since the steam needed for stripping consumes 25–35 % of a power plant's gross output when capture is added.

The current state of the art is **3.0–3.7 GJ/t CO₂** specific reboiler duty for MEA at 90 % capture rate. Aspirational targets with second-generation solvents (KS-1, advanced MEA blends, MDEA-PZ) are 2.0–2.5 GJ/t CO₂.

---

## 2. Reaction chemistry

### Primary and secondary amines (MEA, DEA): zwitterion mechanism

For amines with at least one N-H bond, the dominant absorption pathway is the **zwitterion mechanism** first proposed by Caplow (1968) and adapted by Danckwerts (1979):

**Step 1** — CO₂ attacks the amine nitrogen, forming an unstable zwitterion intermediate:

$$\text{CO}_2 + \text{R}_1\text{R}_2\text{NH} \xrightleftharpoons[k_{-1}]{k_1} \text{R}_1\text{R}_2\text{N}^+\text{HCOO}^-$$

**Step 2** — A nearby base (B = amine, hydroxide, or water) deprotonates the zwitterion:

$$\text{R}_1\text{R}_2\text{N}^+\text{HCOO}^- + \text{B} \xrightleftharpoons[k_{-b}]{k_b} \text{R}_1\text{R}_2\text{NCOO}^- + \text{BH}^+$$

The carbamate ion $\text{R}_1\text{R}_2\text{NCOO}^-$ is the stable end product.

### Overall reaction (MEA)

For MEA (the dominant industrial primary amine), the amine itself is the most abundant base. The two steps combine into:

$$\text{CO}_2 + 2\,\text{MEA} \rightleftharpoons \text{MEACOO}^- + \text{MEAH}^+ \qquad \Delta H = -84.7\ \text{kJ/mol CO}_2$$

The **2:1 stoichiometry** is fundamental: each mole of CO₂ consumes two moles of MEA — one as the nucleophile, one as the base. This sets the maximum theoretical loading at **α = 0.5** (mol CO₂ per mol amine).

In practice, lean loading runs at α = 0.18–0.25 (cannot be lower because the reboiler can't strip the last bit without excessive steam) and rich loading at α = 0.45–0.50 (cannot be higher because the kinetic driving force vanishes near saturation). The **cyclic capacity** is α_rich − α_lean = 0.25–0.30.

### Tertiary amines (MDEA): bicarbonate mechanism

Tertiary amines have no N-H bond and cannot form a zwitterion. Instead, MDEA acts only as a **proton acceptor** for the much slower water-catalyzed CO₂ hydration:

$$\text{CO}_2 + \text{H}_2\text{O} + \text{MDEA} \rightleftharpoons \text{HCO}_3^- + \text{MDEAH}^+ \qquad \Delta H = -55\ \text{kJ/mol CO}_2$$

This has two consequences:
1. **Much lower ΔH** → less reboiler steam (good)
2. **Much slower kinetics** → larger absorber (bad)

Industrial use of MDEA is therefore restricted to applications where the gas residence time is large (selective H₂S removal in refineries) or where MDEA is blended with a fast promoter such as **piperazine** (PZ) — the basis of all "advanced" solvent formulations including KS-1, BASF OASE blue, and Cansolv DC-103.

### Rate law (Aboudheir 2003 termolecular)

Both Caplow's zwitterion and Crooks' termolecular formulations collapse, in the high-amine-concentration regime (industrial absorber conditions), to a second-order rate law:

$$r_{\text{abs}} = k_2 \cdot [\text{CO}_2]_{\text{liq}} \cdot [\text{amine}]_{\text{free}} \quad [\text{mol/(m}^3\cdot\text{s)}]$$

Aboudheir's regression of 196 stopped-flow experiments on MEA at 40–80 °C gives:

$$k_2 = 4.61 \times 10^9 \exp\!\left(-\frac{4412}{T}\right) \quad [\text{m}^3/(\text{kmol}\cdot\text{s})]$$

At 40 °C this is k₂ ≈ 7,000 m³/(kmol·s) — about three orders of magnitude faster than the CO₂/H₂O reaction alone. At 80 °C it rises to 21,000 — which is why pilot plants often run absorbers cooler than the optimal absorber temperature might suggest, since you want to maximize α_rich (favored by lower T) more than k₂.

---

## 3. Mass transfer with reaction

### Two-film theory

At a gas-liquid interface, both films offer mass-transfer resistance:

```
    GAS BULK         |  GAS FILM  |   LIQUID FILM     |    LIQUID BULK
                     |            |                   |
   P_CO2,bulk  →→→→  | P_CO2,int  | C_CO2,int → react |  C_CO2 → ≈ 0
                     |            | with amine        |
                     k_G          | k_L, E            |
```

The overall flux is

$$N_{\text{CO}_2} = K_G \cdot (P_{\text{CO}_2,\text{bulk}} - P_{\text{CO}_2,\text{eq}})$$

with the overall coefficient given by series combination of gas-film and (enhanced) liquid-film resistances:

$$\frac{1}{K_G} = \frac{1}{k_G} + \frac{H_{\text{CO}_2}}{E \cdot k_L}$$

where $H_{\text{CO}_2}$ is the Henry's law constant (CO₂ concentration in liquid / partial pressure) and $E$ is the **enhancement factor**.

### Hatta number — controlling regime

The Hatta number measures whether reaction is fast enough to occur entirely *within* the diffusion film:

$$Ha = \frac{\sqrt{k_2 \cdot C_{\text{amine,free}} \cdot D_{\text{CO}_2}}}{k_L}$$

* **Ha < 0.3**: reaction is slow compared to mass transfer → equilibrium-stage behavior, no enhancement ($E \approx 1$)
* **Ha > 3**: reaction is fast → film reactions enhance transfer ($E > 1$)
* **Ha > 10**: pseudo-first-order regime, reaction goes to completion in the film

For 30 wt% MEA in industrial absorbers: $C_{\text{amine,free}} \approx 4{,}000$ mol/m³, $D_{\text{CO}_2} \approx 1.5\times10^{-9}$ m²/s, $k_L \approx 10^{-4}$ m/s, $k_2(40°C) \approx 7{,}000$ m³/(kmol·s). This gives:

$$Ha = \sqrt{7{,}000 \cdot 4{,}000 \cdot 1.5\times10^{-9}} / 10^{-4} \approx 65$$

→ deep in the fast-reaction enhanced-transfer regime.

### Enhancement factor (DeCoursey 1974)

DeCoursey's algebraic approximation gives $E$ explicitly without needing iterative solution of the Higbie penetration model:

$$E = -\frac{Ha^2}{2(E_\infty - 1)} + \sqrt{\left[\frac{Ha^2}{2(E_\infty - 1)}\right]^2 + \frac{Ha^2 \cdot E_\infty}{E_\infty - 1} + 1}$$

where $E_\infty$ is the **instantaneous-reaction limit** (when reaction is so fast that only the diffusion of fresh amine into the film matters):

$$E_\infty = 1 + \frac{D_{\text{amine}} \cdot C_{\text{amine,free}}}{z \cdot D_{\text{CO}_2} \cdot C_{\text{CO}_2,\text{int}}}$$

with $z$ the stoichiometric ratio ($z=2$ for MEA, $z=1$ for tertiary amines). For 30 wt% MEA at α=0.3, $E_\infty$ is typically 10–50.

The behavior is:
* $Ha \ll \sqrt{E_\infty}$: $E \approx Ha / \tanh(Ha)$ (DanckwertsLaplace solution)
* $Ha \gg \sqrt{E_\infty}$: $E \approx E_\infty$ (instantaneous limit)

For the default industrial case the absorber runs at $E \approx 8\text{-}15$ → reaction provides a ~10× speed-up over physical absorption.

---

## 4. Vapor-liquid equilibrium

### Why VLE matters

VLE sets the **maximum** possible loading at a given gas composition (absorber driving force shrinks to zero at α_eq) and the **minimum** possible loading at the regenerator (stripper bottom equilibrium fixes α_lean).

The function we need is $P_{\text{CO}_2,\text{eq}}(\alpha, T)$ — the equilibrium CO₂ partial pressure over a loaded amine solution.

### Implementation: 4-parameter Posey-style fit

The current model uses a simple closed-form correlation fitted to published 30 wt% MEA data (Jou et al. 1995, Aronu et al. 2011, Hilliard 2008):

$$\ln\left(\frac{P_{\text{CO}_2}}{\text{Pa}}\right) = A + \frac{B}{T} + C \cdot \alpha + \frac{D \cdot \alpha^2}{T}$$

with A = 30.27, B = −8254, C = −2.27, D = 9013.

Sanity-check values (30 wt% MEA):

| α   | T = 40 °C  | T = 80 °C   | T = 120 °C   |
| :-: | :--------: | :---------: | :----------: |
| 0.2 |  ~0.1 kPa  | ~2 kPa      | ~20 kPa      |
| 0.4 |  ~2 kPa    | ~30 kPa     | ~200 kPa     |
| 0.5 |  ~20 kPa   | ~200 kPa    | ~1500 kPa    |

This 4-parameter form is accurate to ±20 % over the operating envelope, which is good enough for first-pass design. For higher accuracy, replace with a full electrolyte-NRTL model (Aronu 2011, Hessen 2013) that accounts for ionic-strength effects.

### Other amines

For DEA, MDEA, PZ the model scales the MEA correlation by a van't Hoff factor based on the ΔH_abs difference:

$$P_{\text{CO}_2,\text{amine}}(α, T) = P_{\text{CO}_2,\text{MEA}}(α, T) \cdot \exp\!\left[\frac{\Delta H_{\text{abs,MEA}} - \Delta H_{\text{abs,amine}}}{R} \left(\frac{1}{T} - \frac{1}{T_{\text{ref}}}\right)\right]$$

Plus an empirical pre-factor (MDEA: ×5, PZ: ×0.4) reflecting the different stoichiometry and binding strength.

---

## 5. Energy balance and reboiler duty

### Why reboiler duty dominates the economics

For a typical 600 MW coal plant with 90 % CO₂ capture:
- Captured CO₂: 400 t/h
- Reboiler duty at 3.5 GJ/t: **1,400 MW thermal** of low-pressure steam
- This is ~30 % of the boiler's gross steam output

Reducing the specific reboiler duty by 0.5 GJ/t corresponds to ~200 MW of saved steam — equivalent to ~60 MW of bottoming-cycle electrical power. The economic difference between a "good" amine system at 3.0 GJ/t and a "poor" one at 4.0 GJ/t is enormous over a 25-year plant life.

### Decomposition of reboiler duty

The reboiler heat must accomplish three things:

#### 1. Drive the desorption reaction (reverse of absorption)

$$Q_{\text{rxn}} = (-\Delta H_{\text{abs}}) \cdot \dot{n}_{\text{CO}_2,\text{cycled}}$$

For MEA this is 84.7 kJ/mol CO₂ → 1.92 GJ/t CO₂. This is the **theoretical minimum** for any amine with the same ΔH; reducing it requires switching to a different solvent class (e.g., tertiary amines at 55 kJ/mol → 1.25 GJ/t theoretical minimum).

#### 2. Generate stripping steam (net water leaving the top of the column)

The stripping action requires water vapor to rise through the column. Most of this water condenses internally on the cold descending rich amine — but the **net water leaving the top** of the column is what the reboiler must vaporize:

$$Q_{\text{steam}} = \Delta H_{\text{vap,H}_2\text{O}} \cdot \dot{n}_{\text{H}_2\text{O,top}}$$

The amount of water leaving the top is set by the top-stage vapor composition:

$$\dot{n}_{\text{H}_2\text{O,top}} = \dot{n}_{\text{CO}_2,\text{cycled}} \cdot \frac{1 - y_{\text{CO}_2,\text{top}}}{y_{\text{CO}_2,\text{top}}}$$

For a well-designed stripper, $y_{\text{CO}_2,\text{top}} \approx 0.6\text{-}0.8$ → steam-to-CO₂ ratio of 0.25–0.7 mol/mol → 0.4–1.0 GJ/t CO₂.

#### 3. Heat the rich amine to T_reb (sensible heat)

After the cross-HX preheats the rich amine to ~110 °C, the reboiler still needs to bring it to 120 °C:

$$Q_{\text{sens}} = \dot{m}_{\text{solv}} \cdot c_p \cdot (T_{\text{reb}} - T_{\text{rich,after-HX}})$$

This is typically 0.5–1.0 GJ/t CO₂ and is **strongly dependent on cross-HX approach temperature**. A 5 K approach gives Q_sens ≈ 0.4 GJ/t; a 15 K approach gives ~1.2 GJ/t. This is why amine plants invest heavily in large cross-HX surface area.

### Total breakdown (default case)

For the 30 wt% MEA default-case simulation at 90 % capture, L/G = 4.0:

| Component  | Value (MW) | Value (GJ/t) | Fraction |
| :--------- | :--------: | :----------: | :------: |
| Reaction   | 235.7      | 1.85         | 54 %     |
| Steam      |  91.2      | 0.71         | 21 %     |
| Sensible   | 108.6      | 0.85         | 25 %     |
| **Total**  | **435.5**  | **3.40**     | 100 %    |

These are within the industrial benchmark (50/30/20 % typical, with total 3.5-3.7 GJ/t).

---

## 6. Packed-column hydrodynamics

### Onda 1968 correlations

For random packings (Pall rings, IMTP, Berl saddles) the canonical correlations are:

**Wetted area** (Onda Eq. 5):

$$\frac{a_w}{a_p} = 1 - \exp\!\left[-1.45 \left(\frac{\sigma_c}{\sigma_L}\right)^{0.75} \!\text{Re}_L^{0.1} \,\text{Fr}_L^{-0.05} \,\text{We}_L^{0.2}\right]$$

with the dimensionless groups defined for the liquid stream. Typical industrial values: $a_w/a_p \approx 0.4\text{-}0.7$.

**Liquid-side coefficient** (Onda Eq. 7):

$$k_L \left(\frac{\rho_L}{\mu_L g}\right)^{1/3} = 0.0051 \,\text{Re}_L^{2/3} \,\text{Sc}_L^{-1/2} \,(a_p d_p)^{0.4}$$

Typical: $k_L \approx 0.5\text{-}2 \times 10^{-4}$ m/s.

**Gas-side coefficient** (Onda Eq. 6):

$$\frac{k_G R T}{a_p D_G} = 5.23 \,\text{Re}_G^{0.7} \,\text{Sc}_G^{1/3} \,(a_p d_p)^{-2}$$

Typical: $k_G \approx 1\text{-}5 \times 10^{-3}$ mol/(m²·s·Pa).

### Packing selection

The current model supports six packings via `PACKING_DB`:

| Packing          | a_p (m²/m³) | ε (–) | d_p (m) | Use                          |
| :--------------- | :---------: | :---: | :-----: | :--------------------------- |
| IMTP #25         | 226         | 0.974 | 0.025   | High a_p, moderate ΔP        |
| **IMTP #50**     | **102**     | **0.978** | **0.050** | **Industry default**       |
| Pall ring 25 mm  | 207         | 0.940 | 0.025   | Older standard               |
| Pall ring 50 mm  | 112         | 0.952 | 0.050   | Older standard               |
| Mellapak 250.Y   | 250         | 0.970 | 0.020   | Structured, low ΔP, high a_p |
| Mellapak 350.Y   | 350         | 0.960 | 0.015   | Most aggressive, modern      |

For amine CO₂ capture, the industry standard is IMTP #50 in the absorber (gives good mass transfer at low pressure drop, easy to clean) and IMTP #25 in the stripper (smaller-scale column, higher a_p OK because pressure drop is less critical).

---

## 7. Reading list (in order of priority)

1. **Kohl & Nielsen, *Gas Purification* 5th ed. (1997)** — the encyclopedic reference. Chapter 2 on amines is the gold standard.
2. **Notz et al., IJGGC 6 (2012) 84–112** — pilot-plant validation data, target dataset for any new amine model.
3. **Aboudheir et al., CES 58 (2003) 5195** — the canonical MEA kinetics paper. Highly readable.
4. **Plaza et al., CEJ 162 (2010) 718** — process modeling and energy integration.
5. **Mac Dowell, Florin, Buchard, ... Energy Environ. Sci. 6 (2013) 2493** — broad review of post-combustion CCS.
6. **IEAGHG 2014/03** — techno-economic baseline.

For the broader context of why CCS matters: IPCC AR6 WG3 Chapter 7 (Energy systems), or the IEA Energy Technology Perspectives 2023.

For a fully-cited deep-dive, see the inline references in the source code and the [`README.md`](../README.md) bibliography.
