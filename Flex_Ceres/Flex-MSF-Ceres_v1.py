"""Ceres Nanoparticle Protocol for Opentrons Flex Robot

Project: Ceres nanoparticle on Flex Robot (V1)
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

# Metadata
metadata = {
    'protocolName': 'Ceres Nanoparticle Protocol for Opentrons Flex Robot (V1)',
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

    # ===== TIMING CONSTANTS =====
    BEAD_BINDING_TIME = 1 if dry_run else 1200  # 20 min
    MAGNET_TIME = 1 if dry_run else 30          # 30 sec Bead pellet time
    WASH_TIME = 1 if dry_run else 120           # 2 min shaking
    ELUTION_TIME = 1 if dry_run else 300        # 5 min elution
    DRYING_TIME = 1 if dry_run else 300         # 5 min drying

    # ===== LIQUID HANDLING CONSTANTS =====
    VOLUME = {
        '': 200,
        'ethanol': 400,
        'wash': 450,
        'supernatant': 500,
        'beads': 30,
        'dnase': 50,
        'elution_water': 30,
    }
    
    ASPIRATE_HEIGHT = 2           # mm above bottom
    AIR_GAP = 20                  # µL
    STANDARD_RATE = 0.6           # Aspiration rate
    
    SHAKE_SPEED = 2000             # RPM
    
    COL_RANGE = range(col_start, col_start + num_col)

    # ===== LABWARE SETUP =====
    
    # Collection plate for eluted RNA
    collection_plate = ctx.load_labware(
        'opentrons_96_wellplate_200ul_pcr_full_skirt', location='D4')
    collected_rna = collection_plate.rows()[0]
    
    # Tip racks
    tips_50_beads = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="A2")
    tips_50_dnase = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="B4")
    tips_50_elution = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_50ul", location="C4")
    tips_1000_solvent = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_1000ul", location="A3")
    tips_1000_supernatant = ctx.load_labware(
        "opentrons_flex_96_filtertiprack_1000ul", location="B3")

    solvent_tips = tips_1000_solvent.rows()[0]
    supernatant_tips = tips_1000_supernatant.rows()[0]
    beads_tips = tips_50_beads.rows()[0]
    dnase_tips = tips_50_dnase.rows()[0]
    elution_tips = tips_50_elution.rows()[0]

    # Pipettes
    p50_multi = ctx.load_instrument(
        "flex_8channel_50", mount="left", tip_racks=[tips_50_beads, tips_50_dnase, tips_50_elution])
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
    #heater_shaker.close_labware_latch()
    
    mag_module = ctx.load_module("magneticBlockV1", location="D1")
    ctx.load_module("thermocyclerModuleV2")  # Loaded but not used in this version
    
    # Reagent and waste plates
    reagent_plate_1 = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="B2")
    waste_plate = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="C2")
    beads_plate = ctx.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", location='D2')
    beads = beads_plate.rows()[0]
    reagent_plate_2 = ctx.load_labware(
        "agilent_6_reservoir_47ml", location="C3")
    
    ctx.load_waste_chute()  # Waste chute available if needed
    
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
        heater_shaker.open_labware_latch()
        ctx.move_labware(sample_plate, mag_module, use_gripper=True)
        ctx.delay(MAGNET_TIME, msg="Pelleting beads on magnetic block.")
    
    def move_to_heater_shaker():
        """Move sample plate back to heater shaker and latch."""
        ctx.move_labware(sample_plate, heater_shaker_adapter, use_gripper=True)
        heater_shaker.close_labware_latch()
    
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
    ctx.comment("Starting RNA purification protocol")

    heater_shaker.open_labware_latch()

    ctx.pause("Please add the sample plate to the heater-shaker module. Press resume when ready.")

    heater_shaker.close_labware_latch()