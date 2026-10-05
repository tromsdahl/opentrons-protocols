# FLEX - Ceres Nanoparticle Isolation Protocol

## Purpose
This purpose of this protocol is to perform the Ceres nanoparticle enrichment on plasma and serum samples for proteomics analysis. The Flex draft protocol is poorly documented, making it difficult to modify or improve. This is a rewrite of the protocol for the Flex robot based off of the MagBead RNA isolation protocol previously generated.

## Development
### 10/05/2026 - Initiation
Initiated writing a Python protocol for Ceres nanoparticle enrichment isolation from plasma/serum samples (assuming 50 uL of starting volume).

#### Deck Layout

|   |  1  |  2  |  3  |  4  |
|---|-----|-----|-----|-----|
| A |  TC |     |     |     |
| B |  TC |     |     |     |
| C |  HS |     |     |     |
| D |  MB |     |     |     |

- TC: Thermocycler
- HS: Heater shaker with deep well adapter
- MP: Magnetic Block

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
1. 