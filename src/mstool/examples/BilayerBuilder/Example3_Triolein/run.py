import mstool
import shutil
shutil.copyfile(mstool.TRIOMAPPING, 'mapping.dat')
shutil.copyfile(mstool.TRIOMARTINI, 'martini.itp')
shutil.copyfile(mstool.TRIOFF,      'ff.xml')

mstool.BilayerBuilder(upper={'POPC':100, 'TRIO':10},
                      lower={'POPC':100, 'TRIO':10},
                      mapping_add='mapping.dat',
                      martini_add='martini.itp',
                      ff_add='ff.xml')

