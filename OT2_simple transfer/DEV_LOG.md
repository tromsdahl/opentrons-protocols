# OT2 - Simple Transfer

## Purpose
This purpose of this protocol is to transfer the peptide digests from a 96-well PCR plate containing the SP3 magnetic beads to a collection plate of the same type. The beads will be pulled to the sides of the wells by the OT2 magnet and the aspiration of the pipette will be offset from the side with the beads and magnet.

We've observed that the OT2 magnet module is able to collect the magnetic beads to a single point much more effectively than the Bravo or Flex magnets that seems to pull to two sides. The Bravo tips are rather thick which causes them to pull up some of the beads during aspiration.

## Development
### 9/30/2026 - Initiation
Initiated writing a Python protocol for transferring 45 uL of peptide digest to a collection plate.

Performed a test for a full plate. It seemed to work well. There was still some residual liquid in the source plate, so it may not be transferring the full 45 uL, but it didn't seem like beads were transferred to the collection plate. This may be due to the height of aspiration being at 2 mm above the bottom (default is 1 mm). But this could be a good thing given that the beads occupy some amount of volume at the bottom of the well. If need be, the protocol could be run twice to be sure all beads have been removed.

### 10/7/2026 - First Use
Used the protocol today. Some changes to make for a version 2.0:
- [ ] Aspirate and dispense a little more slowly so all the volume is transferred
- [ ] Add a touch tip at the end of the dispensing to minimize sample loss and prevent dripping while used tips are being disposed of
- [ ] Reduce the amount of time on the magnet. Five minutes was more than enough time.
