"""Ceres Nanoparticle Protocol for Opentrons Flex Robot

Project: Ceres nanoparticle on Flex Robot (V2)
Date: 10/05/2026
Experiment: Proteomics from Plasma/Serum Samples

Objectives: Perform nanoparticle enrichment using Ceres magnetic isolation protocol
- Nanoparticle enrichment from plasma/serum samples
- Alkylation and reduction of proteins
- Wash and protein aggregation capture
- Trypsin digestion
- TFA quenching and transfer of digested peptides
"""

from opentrons import protocol_api
from opentrons.protocol_api import OFF_DECK

# Metadata
metadata = {
    'protocolName': 'Ceres Nanoparticle Enrichment Protocol (V2)',
    'description': 'Perform nanoparticle enrichment using Ceres magnetic isolation protocol.',
    'author': 'TR'
}

requirements = {'robotType': 'Flex', 'apiLevel': '2.18'}


def add_parameters(parameters):
    """Define protocol parameters for user customization."""
    parameters.add_bool(
        variable_name="dry_run",
        display_name="Dry Run or Sample Run?",
        description="Use shorter incubation times and return tips instead of dropping?",
        default=True
    )

    parameters.add_int(
        variable_name="num_col",
        display_name="Number of columns to run?",
        description="How many columns of a 96 well plate have samples?",
        default=1,
        minimum=1,
        maximum=12
    )

    parameters.add_int(
        variable_name="col_start",
        display_name="Column Start Position",
        description="Column start position on 96 well plate (1-indexed)?",
        default=1,
        minimum=1,
        maximum=12
    )

def run(ctx: protocol_api.ProtocolContext):
    """Execute Nanoparticle purification protocol.
    
    This protocol performs magnetic bead-based nanoparticle purification from plasma/serum samples.
    """
    
    # ===== RUNTIME PARAMETERS =====
    dry_run = ctx.params.dry_run
    num_col = ctx.params.num_col
    col_start = ctx.params.col_start - 1  # Convert to 0-indexed

    # Parameter validation
    if not 1 <= num_col <= 12:
        raise ValueError(f"num_col must be 1-12, got {num_col}")
    if not 0 <= col_start < 12:
        raise ValueError(f"col_start must be 0-11 (internally), got {col_start}")
    if col_start + num_col > 12:
        raise ValueError(f"Sample range exceeds 12 columns: "
                         f"columns {col_start + 1}-{col_start + num_col} invalid")

    # ===== TIMING & TEMP CONSTANTS =====
    MAGNET_TIME = 5 if dry_run else 30          # 30 sec Bead pellet time
    SHAKE_TIME = 5 if dry_run else 30           # 30 sec shaking
    RED_TIME = 5 if dry_run else 1200           # 20 min reduction time
    RED_TEMP = 65                               # 65C Reduction temperature
    RED_COOL_TIME = 5 if dry_run else 300       # 5 min to allow heater-shaker module to cool before alkylation
    ALK_TIME = 5 if dry_run else 1200           # 20 min alkylation time
    DIGEST_TIME = 5 if dry_run else 7200        # 2 hr digestion time
    DIGEST_TEMP = 70                            # 70C Digestion temperature

    # ===== LIQUID HANDLING CONSTANTS =====
    VOLUME = {
        'plasma': 50,
        'nanotrap_buff': 150,
        'bead_A': 66,
        'bead_B': 66,
        'bead_C': 66,
        'H2O_wash': 450,
        'ABC_alk_red': 100,
        'tcep': 5,
        'iodoacetamide': 5,
        'ABC_wash': 450,
        'TEAB': 40,
        'trypsin': 5,
        'TFA_2pct': 5,
    }
    
    ASPIRATE_HEIGHT = 1.5         # mm above bottom
    AIR_GAP = 20                  # µL
    SMALL_AIR_GAP = 5             # µL
    STANDARD_RATE = 0.8           # Aspiration rate
    
    SHAKE_SPEED = 2200            # RPM
    
    COL_RANGE = range(col_start, col_start + num_col)

    # ===== LABWARE SETUP =====
    
    # Collection plate for eluted peptides
    #collection_plate = ctx.load_labware(
    #    'opentrons_96_wellplate_200ul_pcr_full_skirt', location='D2')
    #collected_peptides = collection_plate.rows()[0]
    
    # Tip racks
    tips_50_reduction = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="A2")
    tips_50_alkylation = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="A3")
    tips_50_trypsin = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="B4")
    tips_50_tfa = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="C4")
    tips_50_elution = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="D4")
    tips_1000_solvent = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_1000ul", location="B3")
    tips_1000_supernatant = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_1000ul", location="C3")

    solvent_tips = tips_1000_solvent.rows()[0]
    supernatant_tips = tips_1000_supernatant.rows()[0]
    reduction_tips = tips_50_reduction.rows()[0]
    alkylation_tips = tips_50_alkylation.rows()[0]
    trypsin_tips = tips_50_trypsin.rows()[0]
    tfa_tips = tips_50_tfa.rows()[0]
    elution_tips = tips_50_elution.rows()[0]

    # Pipettes
    p50_multi = ctx.load_instrument(
        "flex_8channel_50", mount="left", tip_racks=[tips_50_reduction, tips_50_alkylation, tips_50_trypsin, tips_50_tfa, tips_50_elution])
    p1000_multi = ctx.load_instrument(
        "flex_8channel_1000", mount="right", 
        tip_racks=[tips_1000_solvent, tips_1000_supernatant])

    # Modules and sample plate
    heater_shaker = ctx.load_module("heaterShakerModuleV1", location="C1")
    heater_shaker_adapter = heater_shaker.load_adapter(
        "opentrons_96_deep_well_adapter")
    sample_plate = heater_shaker_adapter.load_labware(
        "nest_96_wellplate_2ml_deep")
    samples = sample_plate.rows()[0]
    
    mag_module = ctx.load_module("magneticBlockV1", location="D1")
    #ctx.load_module("thermocyclerModuleV2")
    
    # Reagent and waste plates
    reagent_plate_1 = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="B2")
    reagent_plate_2 = ctx.load_labware(
        "axygen_12_reservoir_22000ul", location="C2")
    waste_plate = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="D2")
    
    ctx.load_waste_chute()
    
    # ===== REAGENT CONFIGURATION =====
    REAGENTS = {
        'nanotrap_buff': reagent_plate_2['A1'],
        'beads_A': reagent_plate_2['A2'],
        'beads_B': reagent_plate_2['A3'],
        'beads_C': reagent_plate_2['A4'],
        'ABC': reagent_plate_2['A5'],
        'TCEP': reagent_plate_2['A6'],
        'IAM': reagent_plate_2['A7'],
        'TEAB': reagent_plate_2['A8'],
        'trypsin': reagent_plate_2['A9'],
        'tfa_2pct': reagent_plate_2['A10'],
        'water_wash_1': reagent_plate_1['A1'],
        'water_wash_2': reagent_plate_1['A2'],
        'ABC_wash': reagent_plate_1['A3'],
        'water_wash_3': reagent_plate_1['A4'],
        'water_wash_4': reagent_plate_1['A5'],
    }
    
    WASTE_DISTRIBUTION = {
        'initial_supernatant': waste_plate['A6'],
        'water_wash_1': waste_plate['A1'],
        'water_wash_2': waste_plate['A2'],
        'ABC_wash': waste_plate['A3'],
        'water_wash_3': waste_plate['A4'],
        'water_wash_4': waste_plate['A5'],
        'red_alk_supernatant': reagent_plate_2['A12'],
    }

    # ===== HELPER FUNCTIONS =====
    
    def dispose_tip(pipette):
        """Return or drop tip based on run mode."""
        if dry_run:
            pipette.return_tip()
        else:
            pipette.drop_tip()
    
    def add_reagent(pipette, reagent_well, volume, tip_index=None, dispose_after=False):
        """Add a reagent to all samples in column range.
        
        Args:
            pipette: Pipette instrument to use
            reagent_well: Source well for reagent
            volume: Volume to dispense (µL)
            tip_index: Specific tip index; if None, picks up new tip
        """
        if tip_index is not None:
            pipette.pick_up_tip(solvent_tips[tip_index])
        else:
            pipette.pick_up_tip()
        
        for col in COL_RANGE:
            pipette.aspirate(volume, reagent_well, rate=STANDARD_RATE)
            pipette.aspirate(AIR_GAP, reagent_well.top(10), rate=STANDARD_RATE)
            pipette.dispense(volume + AIR_GAP, samples[col].top(-0.5), rate=STANDARD_RATE)
            pipette.blow_out(samples[col].top(-0.5))
            #pipette.touch_tip(samples[col], radius=0.15, v_offset=0.2, speed=60)
        
        if dispose_after:
            dispose_tip(pipette)
        else:
            pipette.return_tip()
    
    def add_beads(pipette, bead_well, volume, tip_index=None):
        """Add magnetic beads to all samples in column range.
        
        Args:
            pipette: Pipette instrument to use
            bead_well: Source well for beads
            volume: Volume to dispense (µL)
            tip_index: Specific tip index; if None, picks up new tip
        """
        if tip_index is not None:
            pipette.pick_up_tip(solvent_tips[tip_index])
        else:
            pipette.pick_up_tip()
        
        if num_col < 6:
            pipette.mix(3, 125, bead_well)  # Mix beads before aspirating
        elif num_col == 6 or num_col == 7:
            pipette.mix(3, 250, bead_well)  # Mix beads before aspirating
        elif num_col == 8 or num_col == 9:
            pipette.mix(3, 375, bead_well)  # Mix beads before aspirating
        else:
            pipette.mix(3, 500, bead_well)  # Mix beads before aspirating
        
        for col in COL_RANGE:
            pipette.aspirate(volume, bead_well, rate=STANDARD_RATE)
            pipette.aspirate(AIR_GAP, bead_well.top(1), rate=STANDARD_RATE)
            pipette.dispense(volume + AIR_GAP, samples[col].top(-0.5), rate=STANDARD_RATE)
            pipette.blow_out(samples[col].top(-0.5))
            #pipette.touch_tip(samples[col], radius=0.15, v_offset=0.2, speed=60)
        
        dispose_tip(pipette)  # Bead tips are not reused
    
    def remove_supernatant(volume, waste_destinations):
        """Remove supernatant from samples and distribute to waste.
        
        Args:
            volume: Volume to aspirate from each sample (µL)
            waste_destinations: Single well or list of wells for distribution
                If list with 2 elements, alternates between wells
        """
        # Normalize to list format
        if not isinstance(waste_destinations, list):
            waste_destinations = [waste_destinations]
        
        for i, col in enumerate(COL_RANGE):
            pipette = p1000_multi
            pipette.pick_up_tip(supernatant_tips[col])
            pipette.aspirate(volume, samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
            pipette.aspirate(AIR_GAP, samples[col].top(1), rate=STANDARD_RATE)
            
            # Distribute to appropriate waste well
            if len(waste_destinations) == 2:
                waste_well = waste_destinations[i % 2]
            else:
                waste_well = waste_destinations[0]
            
            pipette.dispense(volume + AIR_GAP, waste_well.top(1))
            pipette.blow_out(waste_well)
            pipette.return_tip()  # Supernatant tips are reused
    
    def move_to_magnetic():
        """Move sample plate to magnetic module and pellet beads."""
        heater_shaker.open_labware_latch()
        ctx.move_labware(sample_plate, mag_module, use_gripper=True)
        ctx.delay(MAGNET_TIME, msg="Pelleting beads on magnetic block.")
    
    def move_to_heater_shaker():
        """Move sample plate back to heater shaker and latch."""
        ctx.move_labware(sample_plate, heater_shaker_adapter, use_gripper=True)
        heater_shaker.close_labware_latch()
    
    def wash_and_pellet(reagent_well, volume, waste_dest, msg="", dispose_after=False, tip_index=None):
        """Standard wash cycle: add reagent, shake, move to magnet, pellet, remove supernatant.
        
        Args:
            reagent_well: Source well for wash reagent
            volume: Volume to dispense (µL)
            waste_dest: Destination for supernatant
            msg: Custom delay message
        """
        # Add reagent
        if tip_index is not None:
            p1000_multi.pick_up_tip(solvent_tips[tip_index])
        else:
            p1000_multi.pick_up_tip()
        
        
        for col in COL_RANGE:
            p1000_multi.aspirate(volume, reagent_well, rate=STANDARD_RATE)
            p1000_multi.aspirate(AIR_GAP, reagent_well.top(10), rate=STANDARD_RATE)
            p1000_multi.dispense(volume + AIR_GAP, samples[col].top(1), rate=STANDARD_RATE)
            p1000_multi.blow_out(samples[col].top(1))
            #p1000_multi.touch_tip(samples[col], radius=0.15, v_offset=0.2, speed=60)
        if dispose_after:
            dispose_tip(p1000_multi)
        else:
            p1000_multi.return_tip()
        
        # Shake
        heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
        ctx.delay(SHAKE_TIME, msg=msg or "Mixing beads and reagent.")
        heater_shaker.deactivate_shaker()
        
        # Move and pellet
        move_to_magnetic()
        remove_supernatant(volume, waste_dest)
        move_to_heater_shaker()

    # ===== PROTOCOL START =====
    ctx.comment("Starting Ceres nanoparticle purification protocol")

    heater_shaker.open_labware_latch()

    ctx.pause("Please add the sample plate containing 50 uL of plasma per well to the heater-shaker module. Press resume when ready.")

    heater_shaker.close_labware_latch()

    # add Nanotrap buffer to samples
    add_reagent(p1000_multi, REAGENTS['nanotrap_buff'], VOLUME['nanotrap_buff'], tip_index=0, dispose_after=True)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing Nanotrap buffer with samples.")
    heater_shaker.deactivate_shaker()

    # add beads A, B, C to samples
    add_beads(p1000_multi, REAGENTS['beads_A'], VOLUME['bead_A'], tip_index=1)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing beads A with samples.")
    heater_shaker.deactivate_shaker()

    add_beads(p1000_multi, REAGENTS['beads_B'], VOLUME['bead_B'], tip_index=2)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing beads B with samples.")
    heater_shaker.deactivate_shaker()

    add_beads(p1000_multi, REAGENTS['beads_C'], VOLUME['bead_C'], tip_index=3)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing beads C with samples.")
    heater_shaker.deactivate_shaker()

    # remove supernatant after bead binding
    move_to_magnetic()
    remove_supernatant(VOLUME['plasma'] + VOLUME['nanotrap_buff'] + VOLUME['bead_A'] + VOLUME['bead_B'] + VOLUME['bead_C'], WASTE_DISTRIBUTION['initial_supernatant'])
    move_to_heater_shaker()

    # Wash 1: Water wash
    wash_and_pellet(REAGENTS['water_wash_1'], VOLUME['H2O_wash'], WASTE_DISTRIBUTION['water_wash_1'], msg="Washing beads with water wash 1.", dispose_after=False, tip_index=4)

    # Wash 2: Water wash
    wash_and_pellet(REAGENTS['water_wash_2'], VOLUME['H2O_wash'], WASTE_DISTRIBUTION['water_wash_2'], msg="Washing beads with water wash 2.", dispose_after=False, tip_index=4)

    # Perform reduction
    ## add ABC buffer for reduction
    add_reagent(p1000_multi, REAGENTS['ABC'], VOLUME['ABC_alk_red'], tip_index=5, dispose_after=False)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing ABC buffer with samples.")
    heater_shaker.deactivate_shaker()

    ## add TCEP for reduction
    for col in COL_RANGE:
        p50_multi.pick_up_tip(reduction_tips[col])
        p50_multi.aspirate(VOLUME['tcep'], REAGENTS['TCEP'], rate=STANDARD_RATE)
        p50_multi.aspirate(SMALL_AIR_GAP, REAGENTS['TCEP'].top(10), rate=STANDARD_RATE) # small air gap
        p50_multi.dispense(SMALL_AIR_GAP, samples[col].top(1), rate=STANDARD_RATE) # dispense air gap above sample to prevent bubbles
        p50_multi.dispense(VOLUME['tcep'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.blow_out(samples[col].top(-2))
        #p50_multi.touch_tip(samples[col], radius=0.9, v_offset=-2, speed=60)
        dispose_tip(p50_multi)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing TCEP with samples.")
    heater_shaker.deactivate_shaker()

    heater_shaker.set_and_wait_for_temperature(RED_TEMP)
    ctx.delay(RED_TIME, msg="Reducing proteins with TCEP at 65C.")

    heater_shaker.deactivate_heater()
    ctx.delay(RED_COOL_TIME, msg="Cooling samples to room temperature before alkylation.")

    # Perform alkylation
    ctx.pause("Please add iodoacetamide to the reagent plate. Press resume when ready.")

    ## add IAM for alkylation
    for col in COL_RANGE:
        p50_multi.pick_up_tip(alkylation_tips[col])
        p50_multi.aspirate(VOLUME['iodoacetamide'], REAGENTS['IAM'], rate=STANDARD_RATE)
        p50_multi.aspirate(SMALL_AIR_GAP, REAGENTS['IAM'].top(10), rate=STANDARD_RATE) # small air gap
        p50_multi.dispense(SMALL_AIR_GAP, samples[col].top(1), rate=STANDARD_RATE) # dispense air gap above sample to prevent bubbles
        p50_multi.dispense(VOLUME['iodoacetamide'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.blow_out(samples[col].top(-2))
        #p50_multi.touch_tip(samples[col], radius=0.9, v_offset=-2, speed=60)
        dispose_tip(p50_multi)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing iodoacetamide with samples.")
    heater_shaker.deactivate_shaker()

    ctx.delay(ALK_TIME, msg="Alkylating proteins with iodoacetamide at room temperature.")

    # remove supernatant after reduction and alkylation
    move_to_magnetic()
    remove_supernatant(VOLUME['ABC_alk_red'] + VOLUME['tcep'] + VOLUME['iodoacetamide'], WASTE_DISTRIBUTION['red_alk_supernatant'])
    move_to_heater_shaker()

    # Wash 3: ABC wash
    wash_and_pellet(REAGENTS['ABC_wash'], VOLUME['ABC_wash'], WASTE_DISTRIBUTION['ABC_wash'], msg="Washing beads with ABC buffer.", dispose_after=True, tip_index=5)

    # Wash 4: Water wash
    wash_and_pellet(REAGENTS['water_wash_3'], VOLUME['H2O_wash'], WASTE_DISTRIBUTION['water_wash_3'], msg="Washing beads with water wash 3.", dispose_after=False, tip_index=4)

    # Wash 5: Water wash
    wash_and_pellet(REAGENTS['water_wash_4'], VOLUME['H2O_wash'], WASTE_DISTRIBUTION['water_wash_4'], msg="Washing beads with water wash 4.", dispose_after=True, tip_index=4)

    # Perform trypsin digestion
    ## add TEAB buffer for digestion
    add_reagent(p1000_multi, REAGENTS['TEAB'], VOLUME['TEAB'], tip_index=6, dispose_after=True)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing TEAB buffer with samples.")
    heater_shaker.deactivate_shaker()

    ## exchange p50 tip boxes
    ctx.move_labware(tips_50_reduction, new_location="A4", use_gripper=True)
    ctx.move_labware(tips_50_trypsin, new_location="A2", use_gripper=True)

    ctx.move_labware(waste_plate, new_location=OFF_DECK, use_gripper=False)
    
    ctx.pause("Please add trypsin to the 12-well reagent plate at position A9. Additionally remove the waste plate from D2 and replace with the 96-well collection plate for the eluted peptides. Press resume when ready.")

    # Collection plate for eluted peptides
    collection_plate = ctx.load_labware(
        'opentrons_96_wellplate_200ul_pcr_full_skirt', location='D2')
    collected_peptides = collection_plate.rows()[0]

    ## add trypsin for digestion
    for col in COL_RANGE:
        p50_multi.pick_up_tip(trypsin_tips[col])
        p50_multi.aspirate(VOLUME['trypsin'], REAGENTS['trypsin'], rate=STANDARD_RATE)
        p50_multi.aspirate(SMALL_AIR_GAP, REAGENTS['trypsin'].top(10), rate=STANDARD_RATE) # small air gap
        p50_multi.dispense(SMALL_AIR_GAP, samples[col].top(1), rate=STANDARD_RATE) # dispense air gap above sample to prevent bubbles
        p50_multi.dispense(VOLUME['trypsin'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.blow_out(samples[col].top(-2))
        #p50_multi.touch_tip(samples[col], radius=0.9, v_offset=-2, speed=60)
        dispose_tip(p50_multi)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(5, msg="Mixing trypsin with samples.")
    heater_shaker.deactivate_shaker()

    heater_shaker.set_and_wait_for_temperature(DIGEST_TEMP)
    ctx.delay(DIGEST_TIME, msg="Digesting proteins with trypsin at 70C.")

    heater_shaker.deactivate_heater()

    # Quench digestion with TFA
    ## exchange p50 tip boxes
    ctx.move_labware(tips_50_trypsin, new_location="B4", use_gripper=True)
    ctx.move_labware(tips_50_tfa, new_location="A2", use_gripper=True)
    ctx.move_labware(tips_50_alkylation, new_location="C4", use_gripper=True)
    ctx.move_labware(tips_50_elution, new_location="A3", use_gripper=True)

    ## add TFA to quench digestion
    for col in COL_RANGE:
        p50_multi.pick_up_tip(tfa_tips[col])
        p50_multi.aspirate(VOLUME['TFA_2pct'], REAGENTS['tfa_2pct'], rate=STANDARD_RATE)
        p50_multi.aspirate(SMALL_AIR_GAP, REAGENTS['tfa_2pct'].top(10), rate=STANDARD_RATE) # small air gap
        p50_multi.dispense(SMALL_AIR_GAP, samples[col].top(1), rate=STANDARD_RATE) # dispense air gap above sample to prevent bubbles
        p50_multi.dispense(VOLUME['TFA_2pct'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.blow_out(samples[col].top(-2))
        #p50_multi.touch_tip(samples[col], radius=0.9, v_offset=-2, speed=60)
        dispose_tip(p50_multi)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(SHAKE_TIME, msg="Mixing TFA with samples.")
    heater_shaker.deactivate_shaker()

    # Transfer digested peptides to collection plate
    move_to_magnetic()
    for col in COL_RANGE:
        p50_multi.pick_up_tip(elution_tips[col])
        p50_multi.aspirate(VOLUME['TEAB'] + VOLUME['trypsin'] + VOLUME['TFA_2pct'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.dispense(VOLUME['TEAB'] + VOLUME['trypsin'] + VOLUME['TFA_2pct'], collected_peptides[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.blow_out(collected_peptides[col].top(-5))
        #p50_multi.touch_tip(collected_peptides[col], radius=0.9, v_offset=-2, speed=60)
        dispose_tip(p50_multi)
    
    # Discard supernatant tips
    for col in COL_RANGE:
        p1000_multi.pick_up_tip(supernatant_tips[col])
        dispose_tip(p1000_multi)

    ctx.comment("Ceres nanoparticle purification protocol completed. Please remove the collection plate with digested peptides.")