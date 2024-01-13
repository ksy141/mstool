import mstool
mstool.BilayerBuilder(
        upper={'POPC':20, 'TRIO':2}, 
        lower={'POPC':20, 'TRIO':2}, 
        martini_add='martini.itp',
        mapping_add='mapping.dat',
        ff_add='ff.xml', dt=0.02)
