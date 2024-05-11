from   .universe      import Universe
from   ..lib.distance import distance_matrix, distance_overlap
import numpy as np
import pandas as pd


def solvate(u, out=None, 
            solventdr=4.93, removedr=6.0, waterslab=0.8, waterchain='W', center=True, pbc=False,
            membrane=False):
    
    # if there is a zero, raise an error
    assert np.all(u.dimensions[0:3]), 'check your dimensions'

    Nx  = u.dimensions[0]               // solventdr
    Nx2 = (u.dimensions[0] - waterslab) // solventdr
    if Nx != Nx2: Nx = Nx2

    Ny  = u.dimensions[1]               // solventdr
    Ny2 = (u.dimensions[1] - waterslab) // solventdr
    if Ny != Ny2: Ny = Ny2

    Nz  = u.dimensions[2]               // solventdr
    Nz2 = (u.dimensions[2] - waterslab) // solventdr
    if Nz != Nz2: Nz = Nz2

    if Nx > 3.5 and membrane:
        Nx -= 3
    
    if Ny > 3.5 and membrane:
        Ny -= 3

    xx = np.arange(Nx) * solventdr
    yy = np.arange(Ny) * solventdr
    zz = np.arange(Nz) * solventdr
    
    if center:
        xx -= u.dimensions[0] / 2
        yy -= u.dimensions[1] / 2
        zz -= u.dimensions[2] / 2
    
    # increment in every direction (very useful)
    xyz = np.array(np.meshgrid(xx, yy, zz)).T.reshape(-1, 3)

    data = {'resname': 'W', 'name': 'W', 'resid': np.arange(1, len(xyz) + 1), 
            'chain': waterchain, 
            'x': xyz[:,0], 'y': xyz[:,1], 'z': xyz[:,2]}

    wateru = Universe(data=data)
    wateru.dimensions = u.dimensions
    wateru.cell = u.cell
        
    # waterbox
    if len(u.atoms) == 0:
        if out: wateru.write(out)
        return wateru
    
    pos1 = wateru.atoms[['x','y', 'z']].to_numpy(dtype=np.float64)
    pos2 = u.atoms[['x','y', 'z']].to_numpy(dtype=np.float64)
    dim  = np.asarray(u.dimensions, dtype=np.float64)

    if pbc:
        bA = distance_overlap(pos1, pos2, removedr, u.dimensions)
    else:
        bA = distance_overlap(pos1, pos2, removedr)

    u.atoms = pd.concat([u.atoms, wateru.atoms[bA]], ignore_index=True)

    if out: u.write(out)
    return u


def ionize(u, out=None, qtot=None,
           conc = 0.15, pos='SOD', neg='CLA', waterresname='W', ionchain='',
           posionchain=None, negionchain=None):
    """Add ions at a given concentration.
    conc = ( N_positive_ion / 6.02e23 ) / ( V * 1e-27 ) = (N * 1e4) / (V * 6.02)
    """
    
    # if system is not parameterized, all q = 0;
    if all(u.atoms.type == 'tbd') and not all(u.atoms.resname == waterresname) and not qtot:
        #print("\nYou need to parameterize a system first before adding ions")
        print("Adding ions while assuming the input structure has a net charge of 0")
    
    if qtot:
        qtot = qtot
    else:
        qtot = round(u.atoms.charge.sum())
    
    # instead of using vol, use the number of water * water volume per each bead
    # vol  = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
    vol = len(u.atoms[u.atoms.name == 'W']) * 4.93 ** 3

    if pos in ['MG', 'CAL', 'BAR', 'ZN2', 'CD2']:
        factor = 2
    else:
        factor = 1
    
    # number of positive ions and negative ions
    Npos = int(conc * vol * 6.02 * 1e-4)
    Nneg = qtot + Npos * factor

    # if conc = 0.0 (which will be just adding counterions)
    # Npos = 0 -> neg_add can be negative
    # Nneg should be 0 or positive
    while Nneg < -0.5:
        Npos += 1
        Nneg = qtot + Npos * factor
    
    Nneg = int(Nneg)
    bA = u.atoms.resname == waterresname
    nonwateratoms = u.atoms[~bA]
    wateratoms    = u.atoms[bA].sample(frac=1, ignore_index=True) # a new object

    # Positive ions
    # Note that contrary to usual python slices, both the start and the stop are included
    wateratoms.loc[0:Npos-1, 'charge']  = factor
    wateratoms.loc[0:Npos-1, 'name']    = pos
    wateratoms.loc[0:Npos-1, 'resname'] = pos
    if posionchain:
        wateratoms.loc[0:Npos-1, 'chain']   = posionchain
    else:
        wateratoms.loc[0:Npos-1, 'chain']   = ionchain + '1'
    
    pos_preexisting = u.atoms[u.atoms['name'] == pos]
    if len(pos_preexisting) == 0:
        resid = 0
    else:
        resid = pos_preexisting.resid.max()

    resids = np.arange(resid + 1, resid + 1 + Npos)
    wateratoms.loc[0:Npos-1, 'resid'] = resids


    # Negative ions
    wateratoms.loc[Npos:Npos+Nneg-1, 'charge']  = -1
    wateratoms.loc[Npos:Npos+Nneg-1, 'name']    = neg
    wateratoms.loc[Npos:Npos+Nneg-1, 'resname'] = neg
    if negionchain:
        wateratoms.loc[Npos:Npos+Nneg-1, 'chain']   = negionchain
    else:
        wateratoms.loc[Npos:Npos+Nneg-1, 'chain']   = ionchain + '2'

    neg_preexisting = u.atoms[u.atoms['name'] == neg]
    if len(neg_preexisting) == 0:
        resid = 0
    else:
        resid  = neg_preexisting.resid.max()

    resids = np.arange(resid + 1, resid + 1 + Nneg)
    wateratoms.loc[Npos:Npos+Nneg-1, 'resid'] = resids

    
    # Reset resid of water
    bA2 = wateratoms.resname == waterresname
    resid = len(wateratoms[bA2])
    wateratoms.loc[bA2, 'resid'] = np.arange(1, resid + 1)


    # combine
    u.atoms = pd.concat([nonwateratoms, wateratoms], ignore_index=True)

    # check
    #assert 0 == round(u.atoms.charge.sum()), 'net charge is not zero'
    num_final  = len(u.atoms[u.atoms['name'] == pos])
    conc_final = (num_final * 1e4) / (vol * 6.02)
    print(f"{pos}: {Npos}")
    print(f"{neg}: {Nneg}")
    print(f"{pos}.{neg}: {conc_final:.3f} M\n")
    
    if out: u.write(out)
    return u



def SolvateMartini(structure=None, out=None, t=None,
                   dimensions=None,
                   solventdr=4.93, removedr=6.0, waterslab=0.8, waterchain='W', center=True,
                   conc=0.15, qtot=None, pos='SOD', neg='CLA', waterresname='W', ionchain='ZZ', pbc=True,
                   membrane=False, posionchain=None, negionchain=None):
        
    # make a water box
    if dimensions:
        if isinstance(dimensions, int) or isinstance(dimensions, float):
            dimensions = [dimensions] * 3

        u = Universe()
        u.cell = np.array([[dimensions[0], 0, 0],
                           [0, dimensions[1], 0], 
                           [0, 0, dimensions[2]]])
    
        u.dimensions = np.array([*dimensions[0:3],90,90,90])
    
    # provide a structure
    else:
        if isinstance(structure, Universe):
            u = structure
        else:
            u = Universe(structure)

        if t:
            # make new dimensions / cell based on solute particles
            # exisitng dimensions / cell do not matter
            maxx = u.atoms['x'].max() - u.atoms['x'].min()
            maxy = u.atoms['y'].max() - u.atoms['y'].min()
            maxz = u.atoms['z'].max() - u.atoms['z'].min()
            maxd = max(maxx, maxy, maxz)
            dim  = maxd + 2 * t
            u.dimensions = np.array([dim] * 3 + [90] * 3)
            u.cell = np.array([[dim, 0, 0], [0, dim, 0], [0, 0, dim]])

    solvatedu = solvate(u, solventdr=solventdr, removedr=removedr, waterslab=waterslab,
                        waterchain=waterchain, center=center, pbc=pbc, membrane=membrane)
    if conc == 0.0:
        if out: solvatedu.write(out)
        return solvatedu


    ionizedu  = ionize(solvatedu, out=out, qtot=qtot, conc=conc, pos=pos, neg=neg, 
                       waterresname=waterresname, ionchain=ionchain, posionchain=posionchain, negionchain=negionchain)
    return ionizedu

