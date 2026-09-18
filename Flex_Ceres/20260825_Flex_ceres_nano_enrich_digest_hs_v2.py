"""Ceres nanoparticle protein enrichment and digestion on an Opentrons Flex."""

from math import ceil

from opentrons import protocol_api
from opentrons.types import Point

metadata = {
    "protocolName": "Nanotrap Protein Enrichment and Digestion w/ Heater Shaker",
    "author": "Boren Lin, Opentrons",
}

requirements = {"robotType": "Flex", "apiLevel": "2.25"}


# ===== LIQUID HANDLING CONSTANTS =====
VOLUME = {
    "plasma": 50,
    "enrichment_buffer": 150,
    "nanotrap": 200,
    "water": 500,
    "ammonium_bicarbonate": 100,
    "dtt": 25,
    "iaa": 25,
    "trypsin": 25,
}

# Local aliases keep the liquid-handling sections readable while they are
# migrated to the shared VOLUME mapping.
VOL_PLASMA = VOLUME["plasma"]
VOL_BUFFER_4 = VOLUME["enrichment_buffer"]
VOL_NANOTRAP = VOLUME["nanotrap"]
VOL_H2O = VOLUME["water"]
VOL_ABC = VOLUME["ammonium_bicarbonate"]
VOL_DTT = VOLUME["dtt"]
VOL_IAA = VOLUME["iaa"]
VOL_TRYPSIN = VOLUME["trypsin"]

NUM_COLUMNS = 12
AIR_GAP = 20
ASPIRATE_HEIGHT = 0.2
SOURCE_HEIGHT = 1
WATER_HEIGHT = 2
TIP_COLUMNS_PER_RACK = 12

def add_parameters(parameters):
    """Define runtime parameters for the protocol."""
    parameters.add_bool(
        variable_name="test",
        display_name="Test Run",
        description="Use shortened incubation times and return tips.",
        default=False,
    )
    parameters.add_int(
        variable_name="num_col",
        display_name="Number of Columns",
        description="Number of columns filled with samples (maximum: 12).",
        default=12,
        minimum=1,
        maximum=NUM_COLUMNS,
    )
    parameters.add_int(
        variable_name="pipet_location",
        display_name="P1000 8-ch Position",
        description="How is the P1000 8-channel pipette mounted?",
        default=2,
        choices=[
            {"display_name": "on the right", "value": 1},
            {"display_name": "on the left", "value": 2},
        ],
    )


def run(ctx: protocol_api.ProtocolContext):
    """Run Ceres nanoparticle enrichment followed by protein digestion."""

    global test
    global num_col
    global pipet_location

    global tip_count
    global tip_spare_count

    global open_slot_count


    test = ctx.params.test
    num_col = ctx.params.num_col
    pipet_location = ctx.params.pipet_location

    if not 1 <= num_col <= NUM_COLUMNS:
        raise ValueError(f"num_col must be 1-{NUM_COLUMNS}, got {num_col}")

    #################################

    # One 200 uL tip is used for each plasma and ABC addition, plus three
    # reagent transfers for DTT, IAA, and trypsin.
    num_tipbox = ceil((2 * num_col + 3) / TIP_COLUMNS_PER_RACK)
    
    if pipet_location == 1:
        p1k_8_loc = 'right'
    else:
        p1k_8_loc = 'left'

    tip_count = 0
    tip_spare_count = 0

    open_slot = ['A4', 'B4', 'C4', 'D4']
    open_slot_count = 0


    # deck layout

    hs = ctx.load_module('heaterShakerModuleV1', 'C1')
    hs_adapter = hs.load_adapter('opentrons_96_deep_well_adapter')
    working_plate = hs_adapter.load_labware('nest_96_wellplate_2ml_deep', 'Working Plate') 
    rxn = working_plate.rows()[0][:num_col]

    ctx.load_lid_stack('opentrons_tough_universal_lid', 'B3', 1)

    mag = ctx.load_module('magneticBlockV1', 'D1')


    sample_plate = ctx.load_labware('axygen_96_wellplate_500ul', 'C3', 'Samples')
    sample = sample_plate.rows()[0][:num_col]

    # reagent_12well = ctx.load_labware('nest_12_reservoir_15ml', 'D2', 'Reagents') 
    # buffer = reagent_12well.wells()[:2]
    # nanotrap = reagent_12well.wells()[4:6]

    enrich_re = ctx.load_labware('nest_96_wellplate_2ml_deep', 'C2', 'Reagents for enrichment') 
    digest_re = ctx.load_labware('axygen_96_wellplate_500ul', 'D2', 'Reagents for digestion') 
    
    buffer = enrich_re.rows()[0][:2]
    nanotrap = enrich_re.rows()[0][2:4]
    abc = enrich_re.rows()[0][11]
    
    dtt = digest_re.rows()[0][0]
    iaa = digest_re.rows()[0][1]
    trypsin = digest_re.rows()[0][2]

    water_res = ctx.load_labware('nest_1_reservoir_290ml', 'B2', 'Water') 
    water = water_res.wells()[0]

    waste_res = ctx.load_labware('nest_1_reservoir_290ml', 'B1', 'Liquid Waste') 
    waste = waste_res.wells()[0]

    ctx.load_waste_chute()

    tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')

    if num_tipbox == 4:
        spare_tips = [ctx.load_labware('opentrons_flex_96_tiprack_200ul', slot, '200uL Tips')
                    for slot in ['B4', 'C4', 'D4']]
    elif num_tipbox == 3:
        spare_tips = [ctx.load_labware('opentrons_flex_96_tiprack_200ul', slot, '200uL Tips')
                    for slot in ['B4', 'C4']]
    elif num_tipbox == 2:
        spare_tips = [ctx.load_labware('opentrons_flex_96_tiprack_200ul', slot, '200uL Tips')
                    for slot in ['B4']]

    tips_1k = ctx.load_labware('opentrons_flex_96_tiprack_1000ul', 'A2', '1000uL Tips')
    tips_1k_discard = ctx.load_labware('opentrons_flex_96_tiprack_1000ul', 'A1', '1000uL Tips')

    p1k_8 = ctx.load_instrument('flex_8channel_1000', p1k_8_loc)



    def remove_liquid(vol, drop_tip=0):

        p1k_8.tip_racks = [tips_1k_discard]

        for start in rxn:
            p1k_8.pick_up_tip()
            p1k_8.aspirate(vol+20, start.bottom(z=0.2), rate = 0.1)
            ctx.delay(seconds=2)
            p1k_8.dispense(vol+20, waste.top(z=-2), push_out = 20)
            if drop_tip == 0: p1k_8.return_tip()
            else: p1k_8.drop_tip()

        p1k_8.reset_tipracks()


    def refill_tips(spare):

        global open_slot_count

        ctx.move_labware(tips_200, open_slot[open_slot_count], use_gripper=True)   
        ctx.move_labware(spare_tips[spare], 'A3', use_gripper=True)
        open_slot_count = open_slot_count + 1


    def transfer(vol, start):

        p1k_8.tip_racks = [tips_200]
        p1k_8.pick_up_tip()

        if num_col > 6:

            p1k_8.mix(2, vol*6, start.bottom(z=0.2), rate = 0.5)
            p1k_8.aspirate(vol*6, start.bottom(z=0.2), rate = 0.2)
            ctx.delay(seconds=2)
            p1k_8.air_gap(10)
            p1k_8.dispense(10, rxn[0].top(z=-2))  
            for i in range(6):
                p1k_8.dispense(vol, rxn[i].bottom(z=10), rate = 0.5)  
                ctx.delay(seconds=2)   

            p1k_8.aspirate(vol*(num_col-6), start.bottom(z=0.2), rate = 0.2)
            ctx.delay(seconds=2)
            p1k_8.air_gap(10)
            p1k_8.dispense(10, rxn[6].top(z=-2))  
            for i in range(num_col-6):
                p1k_8.dispense(vol, rxn[6+i].bottom(z=10), rate = 0.5)  
                ctx.delay(seconds=2)  

        else:
            p1k_8.mix(2, vol*num_col, start.bottom(z=0.2), rate = 0.5)
            p1k_8.aspirate(vol*num_col, start.bottom(z=0.2), rate = 0.2)
            ctx.delay(seconds=2)
            p1k_8.air_gap(10)
            p1k_8.dispense(10, rxn[0].top(z=-2))  
            for i in range(num_col):
                p1k_8.dispense(vol, rxn[i].bottom(z=10), rate = 0.5)  
                ctx.delay(seconds=2)             

        if test:
            p1k_8.return_tip() 
        else:
            p1k_8.drop_tip()



    # protocol 



    hs.close_labware_latch()
    ctx.comment('==========Enrichment==========')   

    ### add buffer
    ctx.comment('\n------------Adding Buffer------------\n')

    p1k_8.tip_racks = [tips_1k]
    p1k_8.pick_up_tip()

    if num_col > 6:

        buffer_loc = buffer[0]

        p1k_8.mix(1, VOL_BUFFER_4*6, buffer_loc.bottom(z=1))
        p1k_8.aspirate(VOL_BUFFER_4*6, buffer_loc.bottom(z=1))
        ctx.delay(seconds=1)
        p1k_8.air_gap(20)    
        p1k_8.dispense(20, rxn[0].top(z=-2))           
        for i in range(6):
            p1k_8.dispense(VOL_BUFFER_4, rxn[i].bottom(z=10))
            ctx.delay(seconds=1)

        buffer_loc = buffer[1]

        p1k_8.mix(1, VOL_BUFFER_4*(num_col-6)*0.75, buffer_loc.bottom(z=1))
        p1k_8.aspirate(VOL_BUFFER_4*(num_col-6), buffer_loc.bottom(z=1))
        ctx.delay(seconds=1)
        p1k_8.air_gap(20)    
        p1k_8.dispense(20, rxn[6].top(z=-2))           
        for i in range(num_col-6):
            p1k_8.dispense(VOL_BUFFER_4, rxn[6+i].bottom(z=10))
            ctx.delay(seconds=1)   

    else:

        buffer_loc = buffer[0]
    
        p1k_8.mix(1, VOL_BUFFER_4*num_col*0.75, buffer_loc.bottom(z=1))         
        p1k_8.aspirate(VOL_BUFFER_4*num_col, buffer_loc.bottom(z=1))
        ctx.delay(seconds=1)
        p1k_8.air_gap(20)    
        p1k_8.dispense(20, rxn[0].top(z=-2))           
        for i in range(num_col):
            p1k_8.dispense(VOL_BUFFER_4, rxn[i].bottom(z=10))
            ctx.delay(seconds=1) 

    if test:
        p1k_8.return_tip()   
    else:
        p1k_8.drop_tip()


    ### add plasma
    ctx.comment('\n------------Adding Plasma------------\n')

    p1k_8.tip_racks = [tips_200]
   
    for start, end in zip(sample, rxn):
        p1k_8.pick_up_tip()
        p1k_8.mix(1, VOL_PLASMA, start.bottom(z=1))
        p1k_8.aspirate(VOL_PLASMA, start.bottom(z=1))
        ctx.delay(seconds=1)
        p1k_8.air_gap(10)
        p1k_8.dispense(10, end.top(z=-5))
        p1k_8.mix(2, VOL_PLASMA*0.5, end.bottom(z=1))
        p1k_8.dispense(VOL_PLASMA, end.bottom(z=1), push_out=0)
        ctx.delay(seconds=1)
        p1k_8.blow_out(end.top(z=-5))

        if test: 
            p1k_8.return_tip()
        else:
            p1k_8.drop_tip()

        tip_count = tip_count + 1
        if tip_count == 12:
            refill_tips(tip_spare_count)
            tip_count = 0
            tip_spare_count = tip_spare_count + 1

            del ctx.deck['A3']
            tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')
        


    hs.set_and_wait_for_shake_speed(rpm=2000)
    ctx.delay(seconds=30 if not test else 5)
    hs.deactivate_shaker()


    ### add nanotrap
    ctx.comment('\n------------Adding Nanotrap and Mixing 30 min------------\n')

    p1k_8.tip_racks = [tips_1k]
    p1k_8.pick_up_tip()

    if num_col > 6:

        buffer_loc = nanotrap[0]

        p1k_8.mix(5, VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)         
        p1k_8.aspirate(VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)
        ctx.delay(seconds=1)
        p1k_8.air_gap(20)    
        p1k_8.dispense(20, rxn[0].top(z=-2))           
        for i in range(4):
            p1k_8.dispense(VOL_NANOTRAP, rxn[i].bottom(z=15), rate = 0.2)
            ctx.delay(seconds=1) 
            p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=rxn[i].length/2-0.5)))
            p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=0)))
     
        p1k_8.aspirate(VOL_NANOTRAP*2, buffer_loc.bottom(z=1), rate = 0.2)
        ctx.delay(seconds=1)
        p1k_8.air_gap(20)    
        p1k_8.dispense(20, rxn[4].top(z=-2))           
        for i in range(2):
            p1k_8.dispense(VOL_NANOTRAP, rxn[4+i].bottom(z=15), rate = 0.2)
            ctx.delay(seconds=1) 
            p1k_8.move_to(rxn[4+i].top(z=-2).move(Point(x=rxn[4+i].length/2-0.5)))
            p1k_8.move_to(rxn[4+i].top(z=-2).move(Point(x=0)))


        buffer_loc = nanotrap[1]
        
        if num_col-6 > 4:

            p1k_8.mix(5, VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)         
            p1k_8.aspirate(VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[6].top(z=-2))           
            for i in range(4):
                p1k_8.dispense(VOL_NANOTRAP, rxn[6+i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[6+i].top(z=-2).move(Point(x=rxn[6+i].length/2-0.5)))
                p1k_8.move_to(rxn[6+i].top(z=-2).move(Point(x=0)))        

            p1k_8.aspirate(VOL_NANOTRAP*(num_col-6-4), buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[6+4].top(z=-2))           
            for i in range(num_col-6-4):
                p1k_8.dispense(VOL_NANOTRAP, rxn[6+4+i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[6+4+i].top(z=-2).move(Point(x=rxn[6+4+i].length/2-0.5)))
                p1k_8.move_to(rxn[6+4+i].top(z=-2).move(Point(x=0)))

        else:

            p1k_8.mix(5, VOL_NANOTRAP*(num_col-6)*0.75, buffer_loc.bottom(z=1), rate = 0.2)         
            p1k_8.aspirate(VOL_NANOTRAP*(num_col-6), buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[6].top(z=-2))           
            for i in range(num_col-6):
                p1k_8.dispense(VOL_NANOTRAP, rxn[6+i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[6+i].top(z=-2).move(Point(x=rxn[6+i].length/2-0.5)))
                p1k_8.move_to(rxn[6+i].top(z=-2).move(Point(x=0)))

    else:

        buffer_loc = nanotrap[0]

        if num_col > 4:

            p1k_8.mix(5, VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)         
            p1k_8.aspirate(VOL_NANOTRAP*4, buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[0].top(z=-2))           
            for i in range(4):
                p1k_8.dispense(VOL_NANOTRAP, rxn[i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=rxn[i].length/2-0.5)))
                p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=0)))

            p1k_8.aspirate(VOL_NANOTRAP*(num_col-4), buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[4].top(z=-2))           
            for i in range(num_col-4):
                p1k_8.dispense(VOL_NANOTRAP, rxn[4+i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[4+i].top(z=-2).move(Point(x=rxn[4+i].length/2-0.5)))
                p1k_8.move_to(rxn[4+i].top(z=-2).move(Point(x=0)))

        else:
            p1k_8.mix(5, VOL_NANOTRAP*num_col*0.75, buffer_loc.bottom(z=1), rate = 0.2)         
            p1k_8.aspirate(VOL_NANOTRAP*num_col, buffer_loc.bottom(z=1), rate = 0.2)
            ctx.delay(seconds=1)
            p1k_8.air_gap(20)    
            p1k_8.dispense(20, rxn[0].top(z=-2))           
            for i in range(num_col):
                p1k_8.dispense(VOL_NANOTRAP, rxn[i].bottom(z=15), rate = 0.2)
                ctx.delay(seconds=1) 
                p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=rxn[i].length/2-0.5)))
                p1k_8.move_to(rxn[i].top(z=-2).move(Point(x=0)))

    if test:
        p1k_8.return_tip()
    else:
        p1k_8.drop_tip()


    ### incubate

    for s in [500, 1000, 1500, 2000]:
        hs.set_and_wait_for_shake_speed(rpm=s)
    ctx.delay(seconds=5)
    hs.set_and_wait_for_shake_speed(rpm=1000)
    ctx.delay(seconds=5 if test else 60*30)    
    hs.deactivate_shaker()
    hs.open_labware_latch()


    ### wash
    ctx.comment('\n------------Washing------------\n')

    ctx.move_labware(working_plate, mag, use_gripper=True)
    ctx.delay(seconds=30 if not test else 5)  

    remove_liquid(VOL_BUFFER_4+VOL_PLASMA+VOL_NANOTRAP)  
  
    ctx.move_labware(working_plate, hs_adapter, use_gripper=True)
    hs.close_labware_latch()


    for x in range(1 if test else 2):

        p1k_8.tip_racks = [tips_1k]
        p1k_8.pick_up_tip()

        count = int(num_col//2)
        rest = num_col%2
        
        for i in range(count):

            p1k_8.aspirate(VOL_H2O*2, water.bottom(z=2).move(Point(x=(-1)**i*4.5)))
            ctx.delay(seconds=1)
            for j in range(2):
                p1k_8.dispense(VOL_H2O, rxn[i*2+j].bottom(z=20))
                ctx.delay(seconds=1)  

        if rest != 0:
            p1k_8.aspirate(VOL_H2O, water.bottom(z=2).move(Point(x=4.5)))
            ctx.delay(seconds=1)
            p1k_8.dispense(VOL_H2O, rxn[num_col-1].bottom(z=20))
            ctx.delay(seconds=1) 

        if test:
            p1k_8.return_tip()
        else:
            p1k_8.drop_tip()

        for s in [500, 1000, 1500, 2000]:
            hs.set_and_wait_for_shake_speed(rpm=s)
        ctx.delay(seconds=5)
        hs.set_and_wait_for_shake_speed(rpm=1000)
        ctx.delay(seconds=30 if not test else 5)
        hs.deactivate_shaker()

        hs.open_labware_latch()
        ctx.move_labware(working_plate, mag, use_gripper=True)
        ctx.delay(seconds=30 if not test else 5) 

        remove_liquid(VOL_H2O, 1 if x==1 else 0)  
    
        ctx.move_labware(working_plate, hs_adapter, use_gripper=True)
        hs.close_labware_latch()


    hs.open_labware_latch()
    ctx.pause('Enrichment Complete')
    hs.close_labware_latch()


    ## digestion
    ctx.comment('\n==========Digestion==========\n')

    ### add ammonium bicarbonate
    ctx.comment('\n------------Adding Ammonium Bicarbonate------------\n')

    for end in rxn:

        p1k_8.tip_racks = [tips_200]
        p1k_8.pick_up_tip()

        p1k_8.aspirate(VOL_ABC, abc.bottom(z=0.3), rate = 0.2)
        p1k_8.dispense(VOL_ABC, end.bottom(z=0.2), push_out=0, rate = 0.5)
        p1k_8.mix(5, VOL_ABC*0.75, end.bottom(z=0.2))
        p1k_8.blow_out(end.top(z=-2))
        p1k_8.move_to(end.top(z=-2).move(Point(x=end.length/2-0.5)))

        if test: 
            p1k_8.return_tip()
        else:
            p1k_8.drop_tip()

        tip_count = tip_count + 1
        if tip_count == 12:
            refill_tips(tip_spare_count)
            tip_count = 0
            tip_spare_count = tip_spare_count + 1
    
            del ctx.deck['A3']
            tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')


    ### add reduction reagent and incubate

    transfer(VOL_DTT, dtt) 

    tip_count = tip_count + 1
    if tip_count == 12:
        refill_tips(tip_spare_count)
        tip_count = 0
        tip_spare_count = tip_spare_count + 1

        del ctx.deck['A3']
        tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')
    
    hs.open_labware_latch()
    ctx.move_lid('B3', working_plate, use_gripper=True)
    hs.close_labware_latch()

    hs.set_and_wait_for_temperature(60)  
    for s in [500, 1000, 1500, 2000]:
        hs.set_and_wait_for_shake_speed(rpm=s)
    ctx.delay(seconds=5)

    hs.set_and_wait_for_shake_speed(rpm=1000)
    ctx.delay(minutes=0.1 if test else 60) 
    hs.deactivate_heater()
    ctx.delay(seconds=5 if test else 60*10)
    hs.deactivate_shaker()

    hs.open_labware_latch()
    ctx.move_lid(working_plate, 'B3', use_gripper=True)
    hs.close_labware_latch()


    ### add alkylation reagent and incubate
    ctx.comment('\n------------Adding Alkylation Agent------------\n')

    transfer(VOL_IAA, iaa) 

    tip_count = tip_count + 1
    if tip_count == 12:
        refill_tips(tip_spare_count)
        tip_count = 0
        tip_spare_count = tip_spare_count + 1

        del ctx.deck['A3']
        tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')

    hs.open_labware_latch()
    ctx.move_lid('B3', working_plate, use_gripper=True)
    hs.close_labware_latch()

    for s in [500, 1000, 1500, 2000]:
        hs.set_and_wait_for_shake_speed(rpm=s)
    ctx.delay(seconds=5)

    hs.set_and_wait_for_shake_speed(rpm=1000)
    ctx.delay(minutes=0.1 if test else 30) 
    hs.deactivate_shaker()

    hs.open_labware_latch()
    ctx.move_lid(working_plate, 'B3', use_gripper=True)
    hs.close_labware_latch()


    ### add digestion reagent and incubate
    ctx.comment('\n------------Adding Digestion Reagent------------\n')

    transfer(VOL_TRYPSIN, trypsin) 

    tip_count = tip_count + 1
    if tip_count == 12:
        refill_tips(tip_spare_count)
        tip_count = 0
        tip_spare_count = tip_spare_count + 1

        del ctx.deck['A3']
        tips_200 = ctx.load_labware('opentrons_flex_96_tiprack_200ul', 'A3', '200uL Tips')

    hs.open_labware_latch()
    ctx.move_lid('B3', working_plate, use_gripper=True)
    hs.close_labware_latch()

    hs.set_and_wait_for_temperature(37)  
    for s in [500, 1000, 1500, 2000]:
        hs.set_and_wait_for_shake_speed(rpm=s)
    ctx.delay(seconds=5)

    hs.set_and_wait_for_shake_speed(rpm=1000)
    ctx.delay(minutes=0.1 if test else 60*18) 
    hs.deactivate_heater()
    hs.deactivate_shaker()

