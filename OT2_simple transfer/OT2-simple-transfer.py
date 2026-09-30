"""Simple Transfer Of SP3 Digests from Sample Plate to Collection Plate

Project: Simple Transfer of SP3 Digests (v1.0)
Date: 9/30/2026

Objectives:
- Transfer SP3 digests following SP3 protocol performed on the Agilent Bravo

"""
from opentrons import protocol_api
from opentrons.types import Point

# metadata
metadata = {'protocolName': 'Simple Transfer of SP3 Digests - 96 samples (v1.0)',
            'description': 'Transfer the SP3 digests from the sample plate performed on the Agilent Bravo to a collection plate for full 96-well plate',
            'author': 'TR'}

requirements = {'robotType': 'OT-2', 'apiLevel': '2.15'}

def run(protocol: protocol_api.ProtocolContext):

    # labware definitions
    tips_300 = protocol.load_labware('opentrons_96_tiprack_300ul', 6)
    collection_plate = protocol.load_labware('eppendorftwin.tec_96_wellplate_150ul', 2)

    # pipette definitions
    p300_pipette = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks = [tips_300])
    #p20_pipette = protocol.load_instrument('p20_single_gen2', 'right', tip_racks = [tips_20_1, tips_20_2])

    # module definitions
    mag_mod = protocol.load_module("magnetic module gen2", 4)
    sample_plate = mag_mod.load_labware('eppendorftwin.tec_96_wellplate_150ul')
    samples = sample_plate.rows()[0]

    # sample number definitions
    COL_RANGE = range(0,12) # This is the range of columns that will be processed. For example, if you are processing 12 columns, this will be range(0,12), starting at column 1. If you are processing 8 columns, this will be range(0,8).

    # ==== PROTOCOL START ====
    mag_mod.engage(height_from_base = 5) # Engage the magnetic module at 5 mm height

    protocol.delay(minutes = 5, msg = 'Incubating on the magnetic module for 5 minutes to allow beads to separate.')

    for col in COL_RANGE:
        p300_pipette.pick_up_tip()
        
        if col % 2 == 0:
            offset = -1.5
        else:
            offset = 1.5
        
        source = samples[col].bottom(2).move(Point(x=offset, y=0, z=0))
        source_airgap = samples[col].top(2).move(Point(x=offset, y=0, z=0))
        destination = collection_plate.rows()[0][col].bottom(2)
        destination_blowout = collection_plate.rows()[0][col].top(-2)

        p300_pipette.aspirate(45, source)
        p300_pipette.aspirate(10, source_airgap) # air gap
        p300_pipette.dispense(55, destination)
        p300_pipette.blow_out(destination_blowout)
        p300_pipette.drop_tip()
    
    mag_mod.disengage() # Disengage the magnetic module