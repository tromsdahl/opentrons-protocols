# FLEX - Ceres Nanoparticle Isolation Protocol

## Purpose
This purpose of this protocol is to perform the Ceres nanoparticle enrichment on plasma and serum samples for proteomics analysis. The Flex draft protocol is poorly documented, making it difficult to modify or improve. This is a rewrite of the protocol for the Flex robot based off of the MagBead RNA isolation protocol previously generated.

## Development
### 10/05/2026 - Initiation
Initiated writing a Python protocol for Ceres nanoparticle enrichment isolation from plasma/serum samples (assuming 50 uL of starting volume).

#### Deck Layout

|   |  1  |  2    |  3   |  4  |
|---|-----|-------|------|-----|
| A |  TC |  RT   |  AT  |     |
| B |  TC |  R1   |  SoT | TrT |
| C |  HS |  R2   |  SuT | TT  |
| D |  MB |  CP   |  T   |  ET |

- TC: Thermocycler
- HS: Heater shaker with deep well adapter
- MP: Magnetic Block
- RT: p50 Reduction tips
- R1: 6-well reagent plate
- R2: 12-well reagent plate
- CP: 96-well PCR collection plate
- AT: p50 alkylation tips
- SoT: Solvent tips
- SuT: Supernatant tips
- T: Trash chute
- TrT: p50 Trypsin tips
- TT: p50 TFA tips
- ET: p50 elution tips

#### Manual steps
##### Protein enrichment
1. Add 50 uL plasma
2. Add 150 uL of Nanotrap Buffer 4
3. Add 66 uL of particle A
    1. Vortex for 30 sec
4. Add 66 uL of particle B
    1. Vortex for 30 sec
5. Add 66 uL of particle C
    1. Vortex for 30 sec
6. Vortex for 30 sec to ensure beads are resuspended
7. Incubate for 30 min at RT with shaking
8. Place sample on magnetic rack for 30 sec
9. Remove and discard supernatant
10. Add 500 uL of DI water to pellet
11. Vortex for 30 sec
    1. If needed, centrifuge briefly to pull beads away from sides of tubes or cap
12. Place sample on magnetic rack for 30 sec
13. Remove and discard wash
14. Repeat wash once more (steps 10-13)

##### Protein digestion
1. Add 100 uL 50 mM ABC
2. Vortex 30 sec
3. Add 5 uL of 100 mM TCEP <- *This was changed from the original protocol*
4. Vortex 30 sec
5. Incubate at 65C for 20 min <- *This was changed from the original protocol*
6. Add 5 uL of 187.5 mM iodoacetamide
7. Vortex for 30 sec
8. Incubate at RT for 20 min protected from the light
9. Place on magnetic rack for 30 sec
10. Remove supernatant
11. Wash with 500 uL of ABC
12. Remove wash
13. Wash with 500 uL water 2x, place on magnet and remove wash
14. Add 40 uL of TEAB, vortex.
15. Add 5 uL of trypsin.
16. Incubate at 70C for 2 h.
17. Add 5 uL of 2% TFA
18. Place on magnet and remove digested peptides

### 10/06/2026 - Manual Tests
Tested plasma with buffer and beads. The shaker needs to be set to 2200 rpm to fully resuspend beads. Thirty seconds is sufficient for this. Placing the plate on the magnetic block for 30 sec is also sufficient to pellet particles and form a ring.

Also tested the subsequent washes with water and TEAB where beads start to become clumpy. Once the water is added there is a little bit of a flimy layer at the top where some bead loss may occur but seems unavoidable. Otherwise all the other steps and volumes see to be ok.

Created an initial version 1 of the protocol to test. Note that supernatant tips and supernatant tips are reused until the completion of the protocol.

#### Dry run test
Some things to possibly change:
- [X] The dry run test revealed a missing move to magnetic step to remove the supernatant after the reduction and alkylation. *added*
- [X] Remove pretty much all touch tip steps. They're kind of awkward and clunky. *removed*
- [X] Lower dispensing height to 0 mm or -1 mm so it is within wells (volume is much lower than dispensing height) <- *changed to -0.5 mm*
- [X] Possibly add waste plate where collection plate is and then swap it out to add collection plate
- [X] Add pause step before trypsin digest to add trypsin to reagent well

Made a version 2 that should be ready to test with dry run

### 10/8/2026 - Test run with platelet poor and platelet rich plasma
For some reason there was no digest transferred to the collection plate. It might be that the shaking at 2200 rpm is too much for so little volume (40 uL TEAB + 5 uL trypsin). Or that 2 hrs at 70C would've dried out the solution.

To test:
- [X] 45 uL in wells and shake at 2200 rpm. Is there volume still at the bottom or is it all stuck to the sides?
    - After shaking the plate at 2200 rpm with 45 uL of 50 mM TEAB the solution mostly went back to the bottom. There were a few droplets on the sides but not really enough to explain why there wasn't a final digest solution in the collection plate.

- [X] 45 uL at 70C for 2 hrs uncovered. Does the solution evaporate?
    - Definitely the issue that the temperature is too high with such a low volume. In the test the solution completely evaporated.

---

Some possible improvements:

- Assuming there's solution leftover and it is successfully transferrd to the collection plate. Move the collection plate to the thermocycler at 4C to prevent evaporation and to keep the samples at low temperature.
- Increase volume transfers to be at least 10 uL for TCEP, IAM, Trypsin, and TFA to make the pipetting easier. Possibly increase TEAB volume so the final reaction volume is 100 uL.
- Remove high RPM shake after adding trypsin. Just incubate at room temperature with moderate shaking (beads will have been resuspended after adding TEAB).

