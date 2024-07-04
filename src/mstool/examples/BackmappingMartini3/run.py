import mstool

mstool.Backmap('mix8_419lipids.gro', 
               workdir='mix8_419lipids',
               mapping_add='mapping.dat', ff_add='ff.xml',
               changename_add={':POPX': ':POPC',  ':DLPE': ':DLIPE', ':SAP6': ':SAPI24'},
               turn_off_EMNVT=True)

mstool.Backmap('b2ar_mix8.gro',
               workdir='b2ar_mix8',
               mapping_add='mapping.dat', ff_add='ff.xml',
               changename_add={':POPX': ':POPC',  ':DLPE': ':DLIPE', ':SAP6': ':SAPI24', 
                               ':GLUP': ':GLU',   ':DIPE': ':DLIPE', ':DPSM': ':SSM'},
               turn_off_EMNVT=True)

mstool.Backmap('mix8_3017lipids.gro',
               workdir='mix8_3017lipids',
               mapping_add='mapping.dat', ff_add='ff.xml',
               changename_add={':POPX': ':POPC',  ':DLPE': ':DLIPE', ':SAP6': ':SAPI24'},
               turn_off_EMNVT=True)

mstool.Backmap('mix8_5111lipids.gro',
               workdir='mix8_5111lipids',
               mapping_add='mapping.dat', ff_add='ff.xml',
               changename_add={':POPX': ':POPC',  ':DLPE': ':DLIPE', ':SAP6': ':SAPI24'},
               turn_off_EMNVT=True)

mstool.Backmap('curved.gro',
               workdir='curved',
               mapping_add='mapping.dat', ff_add='ff.xml',
               changename_add={':CDL0': ':TOCL0'},
               turn_off_EMNVT=True)

