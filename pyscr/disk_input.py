import sys
import util

tend = util.yrToSec(sys.argv[1])*1e6
tsave = tend*0.01

file = util.find('diskevol.inp')

util.findAndReplaceLine(file, 'tend', f'tend             =    {util.toSciNot(tend)}\n')
util.findAndReplaceLine(file, 'tsave', f'tsave            =    {util.toSciNot(tsave)}\n')