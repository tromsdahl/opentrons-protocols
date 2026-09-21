"""SP3 Proteomics Protocol for Opentrons Flex Robot

Project: SP3 Proteomics Prep on Flex Robot (V1)
Date: 9/03/2026
Experiment: SP3 on Flex, no manual intervention

Objectives: Perform SP3 proteomics preparation on samples
- Alkylation and reduction
- Wash and protein aggreagtion capture
- Trypsin digestion
- TFA quenching and transfer of digested peptides
"""

from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'SP3 Proteomics Protocol for Opentrons Flex Robot (V1)',
    'description': 'Perform SP3 proteomics preparation on samples.',
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
        description="Number of columns of a 96 well plate that have samples.",
        default=1,
        minimum=1,
        maximum=12
    )

    parameters.add_int(
        variable_name="col_start",
        display_name="Column Start Position",
        description="Column start position on 96 well plate?",
        default=1,
        minimum=1,
        maximum=12
    )

    parameters.add_float(
        variable_name="sample_amt",
        display_name="Starting Sample Volume",
        description="Starting sample volume in µL.",
        default=10,
        minimum=10,
        maximum=50
    )


def run(ctx: protocol_api.ProtocolContext):
    """Execute SP3 proteomics protocol.
    
    This protocol performs SP3 proteomics preparation on samples,
    including alkylation, reduction, wash, protein aggregation capture,
    and trypsin digestion.
    """
    
    # ===== RUNTIME PARAMETERS =====
    dry_run = ctx.params.dry_run
    num_col = ctx.params.num_col
    col_start = ctx.params.col_start - 1  # Convert to 0-indexed
    sample_vol = ctx.params.sample_amt

    # Parameter validation
    if not 1 <= num_col <= 12:
        raise ValueError(f"num_col must be 1-12, got {num_col}")
    if not 0 <= col_start < 12:
        raise ValueError(f"col_start must be 0-11 (internally), got {col_start}")
    if col_start + num_col > 12:
        raise ValueError(f"Sample range exceeds 12 columns: "
                         f"columns {col_start + 1}-{col_start + num_col} invalid")

    # ===== TIMING & TEMP CONSTANTS =====
    MAGNET_TIME = 3 if dry_run else 120         # 2 min Bead pellet time
    WASH_TIME = 3 if dry_run else 120           # 2 min shaking
    RED_ALK_TIME = 3 if dry_run else 60         # 10 min reduction/alkylation time
    DIGESTION_TIME = 3 if dry_run else 7200     # 120 min digestion

    RED_ALK_TEMP = 25 if dry_run else 60        # 60C for reduction/alkylation
    DIGESTION_TEMP = 25 if dry_run else 47      # 47C for digestion

    # ===== LIQUID HANDLING CONSTANTS =====
    VOLUME = {
        'red_alk_buff': ((2 * sample_vol) / 3) + 3 + (1/3),
        'pac_acn' : (2.5 * sample_vol) + 12.5,
        'acn' : 171.5,
        'etoh' : 180,
        'trypsin': 5,
        'teab_digest': 35,
        'teab_resuspend': 100,
        '5pct_tfa' : 10,
    }
    
    ASPIRATE_HEIGHT = 2           # mm above bottom
    AIR_GAP = 20                  # µL
    STANDARD_RATE = 0.6           # Aspiration rate
    
    SHAKE_SPEED = 2000             # RPM
    TRYPSIN_SHAKE_SPEED = 1000     # RPM
    
    COL_RANGE = range(col_start, col_start + num_col)

    # ===== LABWARE SETUP =====
    
    # Collection plate for peptides
    collection_plate = ctx.load_labware(
        'opentrons_96_wellplate_200ul_pcr_full_skirt', location='D2')
    collected_peptides = collection_plate.rows()[0]
    
    # Tip racks
    tips_50_red_alk_buff = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="A2")
    tips_50_pac_acn = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="A3")
    tips_50_trypsin = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="C4")
    tips_200_elution = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_200ul", location="D4")
    tips_200_solvent = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_200ul", location="B3")
    tips_200_supernatant = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_200ul", location="C3")

    red_alk_tips = tips_50_red_alk_buff.rows()[0]
    pac_acn_tips = tips_50_pac_acn.rows()[0]
    trypsin_tips = tips_50_trypsin.rows()[0]
    elution_tips = tips_200_elution.rows()[0]
    solvent_tips = tips_200_solvent.rows()[0]
    supernatant_tips = tips_200_supernatant.rows()[0]
    

    # Pipettes
    p50_multi = ctx.load_instrument(
        "flex_8channel_50", mount="left", tip_racks=[tips_50_red_alk_buff, tips_50_pac_acn, tips_50_trypsin])
    p1000_multi = ctx.load_instrument(
        "flex_8channel_1000", mount="right", 
        tip_racks=[tips_200_solvent, tips_200_supernatant, tips_200_elution])

    # Modules and sample plate
    heater_shaker = ctx.load_module("heaterShakerModuleV1", location="C1")
    heater_shaker_adapter = heater_shaker.load_adapter(
        "opentrons_96_pcr_adapter")
    sample_plate = heater_shaker_adapter.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt")
    samples = sample_plate.rows()[0]
    
    mag_module = ctx.load_module("magneticBlockV1", location="D1")
    thermocycler = ctx.load_module("thermocyclerModuleV2")
    
    # Reagent and waste plates
    reagent_plate_1 = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="B2")
    waste_plate = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="C2")
    beads_plate = ctx.load_labware( # to be removed
        "opentrons_96_wellplate_200ul_pcr_full_skirt", location='D2')
    beads = beads_plate.rows()[0] # to be removed
    reagent_plate_2 = ctx.load_labware( # to be removed
        "agilent_6_reservoir_47ml", location="C3")
    
    ctx.load_waste_chute()
    
    # ===== REAGENT CONFIGURATION =====
    REAGENTS = {
        'lysis': reagent_plate_1['A1'],
        'ethanol': reagent_plate_1['A2'],
        'wash_1': reagent_plate_1['A3'],
        'wash_2': reagent_plate_1['A4'],
        'ethanol_1': reagent_plate_1['A5'],
        'ethanol_2': reagent_plate_1['A6'],
        'dnase': reagent_plate_2['A1'],
        'rna_buffer': reagent_plate_2['A2'],
        'dnase_wash_1': reagent_plate_2['A3'],
        'dnase_wash_2': reagent_plate_2['A4'],
        'elution': reagent_plate_2['A5'],
    }
    
    WASTE_DISTRIBUTION = {
        'initial_supernatant': [waste_plate['A1'], waste_plate['A2']],
        'wash_1': waste_plate['A3'],
        'wash_2': waste_plate['A4'],
        'ethanol_1': waste_plate['A5'],
        'ethanol_2': waste_plate['A6'],
        'dnase_supernatant': [reagent_plate_2['A1'], reagent_plate_2['A2']],
    }

    # ===== HELPER FUNCTIONS =====
    
    def dispose_tip(pipette):
        """Return or drop tip based on run mode."""
        if dry_run:
            pipette.return_tip()
        else:
            pipette.drop_tip()
    
    def add_reagent(pipette, reagent_well, volume, tip_index=None):
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
            pipette.dispense(volume + AIR_GAP, samples[col].top(1), rate=STANDARD_RATE)
            pipette.blow_out()
            pipette.touch_tip(radius=0.2)
        
        dispose_tip(pipette)
    
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
            
            pipette.dispense(volume + AIR_GAP, waste_well)
            pipette.return_tip()  # Supernatant tips are reused
    
    def move_to_magnetic():
        """Move sample plate to magnetic module and pellet beads."""
        ctx.move_labware(sample_plate, mag_module, use_gripper=True)
        ctx.delay(MAGNET_TIME, msg="Pelleting beads on magnetic block.")
    
    def move_to_heater_shaker():
        """Move sample plate back to heater shaker and latch."""
        heater_shaker.open_labware_latch()
        ctx.move_labware(sample_plate, heater_shaker_adapter, use_gripper=True)
        heater_shaker.close_labware_latch()

    def move_to_thermocycler(temp, time):
        """Move sample plate to thermocycler for incubations and use heated lid to prevent evaporation."""
        thermocycler.open_lid()
        ctx.move_labware(sample_plate,thermocycler, use_gripper=True)
        thermocycler.close_lid()
        thermocycler.set_temperature(temp)
        ctx.delay(time, msg=f"Incubating at {temp}°C for {time} seconds.")
    
    def wash_and_pellet(reagent_well, volume, waste_dest, msg=""):
        """Standard wash cycle: add reagent, shake, move to magnet, pellet, remove supernatant.
        
        Args:
            reagent_well: Source well for wash reagent
            volume: Volume to dispense (µL)
            waste_dest: Destination for supernatant
            msg: Custom delay message
        """
        # Add reagent
        p1000_multi.pick_up_tip()
        for col in COL_RANGE:
            p1000_multi.aspirate(volume, reagent_well, rate=STANDARD_RATE)
            p1000_multi.aspirate(AIR_GAP, reagent_well.top(10), rate=STANDARD_RATE)
            p1000_multi.dispense(volume + AIR_GAP, samples[col].top(1), rate=STANDARD_RATE)
        dispose_tip(p1000_multi)
        
        # Shake
        heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
        ctx.delay(WASH_TIME, msg=msg or "Mixing beads and reagent.")
        heater_shaker.deactivate_shaker()
        
        # Move and pellet
        move_to_magnetic()
        remove_supernatant(volume, waste_dest)
        move_to_heater_shaker()

    # ===== PROTOCOL START =====
    ctx.comment("Starting SP3 protocol")

    heater_shaker.open_labware_latch()

    ctx.pause("Please add the sample plate to the heater-shaker module. Press resume when ready.")

    heater_shaker.close_labware_latch()
    
    # Add RNA lysis buffer to samples
    add_reagent(p1000_multi, REAGENTS['lysis'], VOLUME['lysis_buffer'], tip_index=0)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(WASH_TIME, msg="Shaking samples with RNA Lysis Buffer for 2 min.")
    heater_shaker.deactivate_shaker()

    # Add Ethanol
    add_reagent(p1000_multi, REAGENTS['ethanol'], VOLUME['ethanol'], tip_index=1)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(WASH_TIME, msg="Shaking samples with RNA Lysis Buffer & EtOH for 2 min.")
    heater_shaker.deactivate_shaker()

    # Add MagBinding Beads
    for col in COL_RANGE:
        p50_multi.pick_up_tip(beads_tips[col])
        p50_multi.mix(10, 20, beads[col].bottom(ASPIRATE_HEIGHT), rate=0.4)
        p50_multi.aspirate(VOLUME['beads'], beads[col], rate=0.4)
        p50_multi.dispense(VOLUME['beads'], samples[col], rate=0.4)
        dispose_tip(p50_multi)
    
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(BEAD_BINDING_TIME, msg="Mixing beads for binding.")
    heater_shaker.deactivate_shaker()

    # First magnetic separation
    move_to_magnetic()
    remove_supernatant(830, WASTE_DISTRIBUTION['initial_supernatant'])
    move_to_heater_shaker()

    # ===== WASH PHASE =====
    # Wash 1
    wash_and_pellet(REAGENTS['wash_1'], VOLUME['wash'],
                    WASTE_DISTRIBUTION['wash_1'], msg="Mixing beads and Wash 1.")

    # Wash 2
    wash_and_pellet(REAGENTS['wash_2'], VOLUME['wash'],
                    WASTE_DISTRIBUTION['wash_2'], msg="Mixing beads and Wash 2.")

    # Ethanol Wash 1
    wash_and_pellet(REAGENTS['ethanol_1'], VOLUME['wash'],
                    WASTE_DISTRIBUTION['ethanol_1'], msg="Mixing beads and EtOH Wash 1.")

    # Ethanol Wash 2
    wash_and_pellet(REAGENTS['ethanol_2'], VOLUME['wash'],
                    WASTE_DISTRIBUTION['ethanol_2'], msg="Mixing beads and EtOH Wash 2.")

    # Move beads plate to A4 for DNase I treatment if needed
    ctx.move_labware(beads_plate, "A4", use_gripper=True)

    # ===== ELUTION PHASE =====
    
    # Dry beads while still on magnetic block
    move_to_magnetic()
    ctx.delay(DRYING_TIME, msg="Drying beads on magnetic block.")
    move_to_heater_shaker()

    # Add DNase/RNase-Free water to elute RNA
    add_reagent(p1000_multi, REAGENTS['elution'], VOLUME['elution_water'], tip_index=7)
    heater_shaker.set_and_wait_for_shake_speed(SHAKE_SPEED)
    ctx.delay(ELUTION_TIME, msg="Eluting RNA from beads.")
    heater_shaker.deactivate_shaker()

    # Final magnetic separation
    heater_shaker.open_labware_latch()
    ctx.move_labware(sample_plate, mag_module, use_gripper=True)
    heater_shaker.close_labware_latch()
    ctx.delay(MAGNET_TIME, msg="Pelleting beads on magnetic block.")

    # Transfer eluted RNA to collection plate
    # Move collection plate and tip boxes to appropriate positions for elution transfer
    ctx.move_labware(collection_plate, "D2", use_gripper=True)
    ctx.move_labware(tips_50_beads, "D4", use_gripper=True)
    ctx.move_labware(tips_50_elution, "A2", use_gripper=True)

    for col in COL_RANGE:
        p50_multi.pick_up_tip(elution_tips[col])
        p50_multi.aspirate(VOLUME['elution_water'], samples[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        p50_multi.dispense(VOLUME['elution_water'], collected_rna[col].bottom(ASPIRATE_HEIGHT), rate=STANDARD_RATE)
        dispose_tip(p50_multi)
    
    # Discard used supernatant tips to clean up
    for col in COL_RANGE:
        p1000_multi.pick_up_tip(supernatant_tips[col])
        dispose_tip(p1000_multi)
    
    ctx.comment("RNA purification protocol complete. Samples are now in collection plate.")
