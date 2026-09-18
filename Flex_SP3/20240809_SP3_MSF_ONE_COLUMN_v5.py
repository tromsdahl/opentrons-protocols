# Project: SP3 for Proteomics Prep (V5)
# Date: 8/9/2024
# Experiment: MSF SP3 on Opentrons OT2

# load libraries
from opentrons import protocol_api

# metadata
metadata = {'protocolName': 'SP3 for 8 samples (Columns 1) V5',
            'description': 'Perform SP3 proteomics prep for 8 samples for 1 column of a 96-well plate.',
            'author': 'TR'}

requirements = {'robotType': 'OT-2', 'apiLevel': '2.14'}

# run function
def run(protocol: protocol_api.ProtocolContext):

    # labware definitions
    #samples = protocol.load_labware('eppendorftwin.tec_96_wellplate_150ul', 2)
    tips_300_1 = protocol.load_labware('opentrons_96_tiprack_300ul', 9)
    tips_300_2 = protocol.load_labware('opentrons_96_tiprack_300ul', 6)
    tips_20_1 = protocol.load_labware('opentrons_96_tiprack_20ul', 8)
    tips_20_2 = protocol.load_labware('opentrons_96_tiprack_20ul', 11)
    solvents = protocol.load_labware('axygen_12_reservoir_22000ul', 5)
    micro_rack = protocol.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', 1)
    collection_plate = protocol.load_labware('eppendorftwin.tecpcrstrip_96_wellplate_150ul', 2)

    # pipette definitions
    p300_pipette = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks = [tips_300_1, tips_300_2])
    p20_pipette = protocol.load_instrument('p20_single_gen2', 'right', tip_racks = [tips_20_1, tips_20_2])

    # module definitions
    mag_mod = protocol.load_module("magnetic module gen2", 4)
    temp_mod = protocol.load_module("temperature module", 7)

    # sample number definitions
    num_columns = 1 # How many columns of samples are you processing?
    column_start = 0 # Where on the plate are the samples starting from (e.g. first column, second, etc.)? Note that the first column is 0, second is 1, etc.

    # set temp module to 60C
    temp_mod.set_temperature(60)

    # REDUCTION AND ALKYLATION
    samples_temp = temp_mod.load_labware('eppendorftwin.tec_96_wellplate_150ul')

    samples = samples_temp.wells()

    # reduction
    for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.aspirate(5, micro_rack['A2'], rate=0.75) # DTT
        #p20_pipette.aspirate(5, solvents['A1'].top(10), rate=0.5) # air gap
        p20_pipette.dispense(5, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        p20_pipette.mix(2, 10, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.5, v_offset=-10, speed=20)
        p20_pipette.drop_tip()
    
    protocol.delay(minutes=30, msg="Incubating for 30 min at 60C for reduction.")

    temp_mod.set_temperature(25)
    
    # alkylation
    for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.aspirate(5, micro_rack['A3'], rate=0.75) # CAA
        #p20_pipette.aspirate(5, solvents['A2'].top(10), rate=0.5) # air gap
        p20_pipette.dispense(5, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        p20_pipette.mix(2, 10, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.5, v_offset=-10, speed=20)
        p20_pipette.drop_tip()


    protocol.delay(minutes=30, msg="Incubating for 30 min at 25C for alkylation.")


    # PROTEIN BINDING
    """# add magnetic beads
    for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.mix(5, 10, micro_rack['A1']) # resuspend the beads
        p20_pipette.touch_tip(radius=0.95, speed=40) # remove any magnetic bead solution from tip
        p20_pipette.aspirate(5, micro_rack['A1'], rate=0.75) # Magnetic beads
        #p20_pipette.aspirate(5, solvents['A2'].top(10), rate=0.5) # air gap
        p20_pipette.dispense(5, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        p20_pipette.mix(2, 10, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.5, v_offset=-9, speed=20)
        p20_pipette.drop_tip()"""

    # add acetonitrile
    
    samples_3 = samples_temp.rows()[0]
    
    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(30, solvents['A3'], rate=0.5) # ACN
        p300_pipette.aspirate(10, solvents['A3'].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, samples_3[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(30, samples_3[i + column_start].bottom(2), rate=0.5)
        p300_pipette.mix(2, 30, samples_3[i + column_start].bottom(3), rate = 0.5)
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.75, v_offset=-7, speed=40)
        p300_pipette.drop_tip()
    
    """for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.aspirate(15, solvents['A3'], rate=0.75) # ACN
        p20_pipette.aspirate(5, solvents['A3'].top(10), rate=0.75) # air gap
        p20_pipette.dispense(5, samples[i + (column_start * 8)].top(1), rate=0.75) # expel air gap before dispensing solution to prevent air bubbles
        p20_pipette.dispense(15, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        p20_pipette.mix(2, 20, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.75, v_offset=-7, speed=40)
        p20_pipette.drop_tip()"""

    protocol.pause("Please move plate from temperature module to magnetic module. Resume when plate has been moved.")

    samples_mag = mag_mod.load_labware('eppendorftwin.tec_96_wellplate_150ul')

    samples_2 = samples_mag.rows()[0]

    mag_mod.engage(height_from_base=5)

    protocol.delay(minutes=5, msg="Incubating on magnet for 5 min")

    # transfer supernatant to waste column of solvent plate
    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(65, samples_2[i + column_start].bottom(2), rate=0.5) # supernatant
        p300_pipette.aspirate(10, samples_2[i + column_start].top(10), rate=0.5) # air gap
        p300_pipette.dispense(75, solvents['A12'].bottom(2), rate=0.5)
        p300_pipette.blow_out()
        p300_pipette.drop_tip()

    # ETHANOL WASH #1
    mag_mod.disengage()

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(200, solvents['A4'], rate=0.5) # EtOH
        p300_pipette.aspirate(10, solvents['A4'].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, samples_2[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(200, samples_2[i + column_start].bottom(2), rate=0.5)
        p300_pipette.mix(10, 100, samples_2[i + column_start].bottom(6), rate = 0.5)
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.95, speed=40)
        p300_pipette.drop_tip()
    
    mag_mod.engage(height_from_base=5)

    protocol.delay(minutes=5, msg="Incubating on magnet for 5 min")

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(230, samples_2[i + column_start].bottom(2), rate=0.5) # supernatant
        p300_pipette.aspirate(10, samples_2[i + column_start].top(10), rate=0.5) # air gap
        p300_pipette.dispense(240, solvents['A12'].bottom(2), rate=0.5)
        p300_pipette.blow_out()
        p300_pipette.drop_tip()

    
    # ETHANOL WASH #2
    mag_mod.disengage()

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(200, solvents['A4'], rate=0.5) # EtOH
        p300_pipette.aspirate(10, solvents['A4'].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, samples_2[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(200, samples_2[i + column_start].bottom(2), rate=0.5)
        p300_pipette.mix(10, 100, samples_2[i + column_start].bottom(6), rate = 0.5)
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.95, speed=40)
        p300_pipette.drop_tip()
    
    mag_mod.engage(height_from_base=5)

    protocol.delay(minutes=5, msg="Incubating on magnet for 5 min")

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(230, samples_2[i + column_start].bottom(2), rate=0.5) # supernatant
        p300_pipette.aspirate(10, samples_2[i + column_start].top(10), rate=0.5) # air gap
        p300_pipette.dispense(240, solvents['A12'].bottom(2), rate=0.5)
        p300_pipette.blow_out()
        p300_pipette.drop_tip()

    mag_mod.disengage()


    # ACN WASH
    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(171.5, solvents['A3'], rate=0.5) # ACN
        p300_pipette.aspirate(10, solvents['A3'].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, samples_2[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(171.5, samples_2[i + column_start].bottom(2), rate=0.5)
        p300_pipette.mix(10, 100, samples_2[i + column_start].bottom(3), rate = 0.5)
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.95, speed=40)
        p300_pipette.drop_tip()
    
    mag_mod.engage(height_from_base=5)

    protocol.delay(minutes=5, msg="Incubating on magnet for 5 min")

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(230, samples_2[i + column_start].bottom(2), rate=0.5) # supernatant
        p300_pipette.aspirate(10, samples_2[i + column_start].top(10), rate=0.5) # air gap
        p300_pipette.dispense(240, solvents['A12'].bottom(2), rate=0.5)
        p300_pipette.blow_out()
        p300_pipette.drop_tip()

    mag_mod.disengage()
    
    protocol.pause("Place plate on temperature module.")

    # ON BEAD DIGESTION
    # add ammonium bicarbonate
    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(35, solvents['A5'], rate=0.5) # ABC
        p300_pipette.aspirate(10, solvents['A5'].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, samples_3[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(35, samples_3[i + column_start].bottom(2), rate=0.5)
        p300_pipette.mix(2, 20, samples_3[i + column_start].bottom(3), rate = 0.5)
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.75, v_offset=-7, speed=40)
        p300_pipette.drop_tip()

    # add trypsin
    for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.aspirate(5, micro_rack['A4'], rate=0.75) # trypsin
        #p20_pipette.aspirate(5, solvents['A3'].top(10), rate=0.5) # air gap
        p20_pipette.dispense(5, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        p20_pipette.mix(2, 10, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.5, v_offset=-8, speed=20)
        p20_pipette.drop_tip()
    
    temp_mod.set_temperature(47)

    #protocol.pause("Please shake the plate on thermomixer for 30 sec at room temp to resuspend beads, then return plate to temperature module that is set to 47C.")

    protocol.delay(minutes=120, msg="Incubate plate for 2 hr at 47C.")

    # ACIDIFICATION WITH 10% TFA
    for i in range(8*num_columns):
        p20_pipette.pick_up_tip()
        p20_pipette.aspirate(10, micro_rack['A5'], rate=0.75) # 10% TFA
        #p20_pipette.aspirate(5, solvents['A3'].top(10), rate=0.5) # air gap
        p20_pipette.dispense(10, samples[i + (column_start * 8)].bottom(2), rate=0.75)
        #p20_pipette.mix(5, 20, samples[i + (column_start * 8)].bottom(3), rate = 0.75)
        p20_pipette.blow_out()
        p20_pipette.touch_tip(radius=0.5, v_offset=-8, speed=20)
        p20_pipette.drop_tip()

    temp_mod.set_temperature(25)

    # Add additional TEAB to account for loss due to evaporation
    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(100, solvents['A5'], rate = 0.5) # add TEAB
        p300_pipette.aspirate(10, solvents['A5'].top(10), rate = 0.5) # air gap
        p300_pipette.dispense(10, samples_3[i + column_start].top(1), rate = 0.5) # expel air gap
        p300_pipette.dispense(100, samples_3[i + column_start].bottom(2), rate = 0.5)
        p300_pipette.mix(5, 80, samples_3[i + column_start].bottom(3)) # mix with full pipetting speed
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.75, v_offset=-7, speed=40)
        p300_pipette.drop_tip()

    # TRANSFER TO COLLECTION PLATE
    protocol.pause("Remove plate from temperature module and place back on to the magnetic module.")

    mag_mod.engage(height_from_base = 5)

    protocol.delay(minutes=5, msg="Incubating on magnet for 5 min")

    collect = collection_plate.rows()[0]

    for i in range(num_columns):
        p300_pipette.pick_up_tip()
        p300_pipette.aspirate(160, samples_2[i + column_start].bottom(2), rate=0.5) # digested samples
        p300_pipette.aspirate(10, samples_2[i + column_start].top(10), rate=0.5) # air gap
        p300_pipette.dispense(10, collect[i + column_start].top(1), rate=0.5) # expel air gap before dispensing solution to prevent air bubbles
        p300_pipette.dispense(160, collect[i + column_start].bottom(1), rate=0.5)
        #p300_pipette.mix(2, 30, collect[i + column_start].bottom(1), rate = 0.5) # mix the digested peptides with TFA to inactivate trypsin
        p300_pipette.blow_out()
        p300_pipette.touch_tip(radius=0.75, v_offset=-6, speed=20)
        p300_pipette.drop_tip()
    
    # disengage magnet module and deactivate the temperature module
    mag_mod.disengage()
    temp_mod.deactivate()