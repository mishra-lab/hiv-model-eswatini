import sys
import model.scenario
from model.scenario import calibrate

model.scenario.uid = '2022-04-20'
# model.scenario.N['cal'] = 1000 # DEBUG
model.scenario.N['b'] = int(sys.argv[1])

calibrate.run()
# calibrate.merge()
# calibrate.rerun()

