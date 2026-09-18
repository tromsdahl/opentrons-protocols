# Project: SP3 for Proteomics Prep on Flex Robot (V1)
# Date: 6/16/2026
# Experiment: MSF SP3 on Opentrons OT2

# Objectives: Perform SP3 protocol on the Opentrons Flex robot without manual intervention (i.e. using the gripper to move sample plate on and off of magnet plate)

# load libraries
from opentrons import protocol_api

# metadata
metadata = {'protocolName': 'SP3 for 8 samples (Column 1) V1',
            'description': 'Perform SP3 proteomics prep for 8 samples for 1 column of a 96-well plate.',
            'author': 'TR'}

requirements = {'robotType': 'Flex', 'apiLevel': '2.15'}

# parameters to modify the protocol
def add_parameters(parameters):
    # reduces incubation times to a few seconds and replaces tips rather than disposes them
    parameters.add_bool(
        variable_name = "dry_run",
        display_name = "Dry Run or Sample Run?",
        description = "Do you want to run a dry run?",
        default = True
    )

    # number of columns with samples
    parameters.add_int(
        variable_name = "num_col",
        display_name = "Number of columnes to run?",
        description = "How many columns of a 96 well plate have samples?",
        default = 1,
        minimum = 1,
        maximum = 12
    )

    # column starting position of 96 well plate
    parameters.add_int(
        variable_name = "col_start",
        display_name = "Column Start Position",
        description = "Column start position on 96 well plate?",
        default = 1,
        minimum = 1,
        maximum = 12
    )

    # variable starting volume of sample (minimum = 10 uL, maximum = )
    parameters.add_float(
        variable_name = "sample_vol",
        display_name = "Sample Volume",
        description = "How much sample is used as the starting volume?",
        default = 10,
        minimum = 5,
        maximum = 45
    )


# run function
def run(protocol: protocol_api.ProtocolContext):

    # runtime parameters
    dry_run = protocol.params.dry_run # Boolean: True or False
    num_col = protocol.params.num_col # integer
    col_start = protocol.params.col_start - 1 # integer (subtract 1 to 0-index for Python slicing)
    sample_vol = protocol.params.sample_vol # float

    # labware definitions

    # collection plate
    collection_plate = protocol.load_labware('eppendorftwin.tecpcrstrip_96_wellplate_150ul', location = 'D2')

    # sample plate
    sample_plate = protocol.load_labware('eppendorftwin.tec_96_wellplate_150ul', location = 'C2')
    
    # tips
    tips_50_1 = protocol.load_labware("opentrons_flex_96_filtertiprack_50ul", location = "")

    # pipette definitions
    p50_multi = protocol.load_instrument("flex_8channel_50", mount = "left")
    p1000_multi = protocol.load_instrument("flex_8channel_1000", mount = "right")
    gripper = protocol.load_instrument("flex_gripper")

    # module definitions
    therm_cyc = protocol.load_module("thermocyclerModuleV2")
    temp_mod = protocol.load_module("temperature module gen2", location = "C1")
    
    temp_mod.set_temperature(4)

    mag_plate = protocol.load_module("magneticBlockV1", location = "D1")
    trash_chute = protocol.load_waste_chute()

    # reagent plate on temp module
    reagent_plate = temp_mod.load_labware('axygen_12_reservoir_22000ul')