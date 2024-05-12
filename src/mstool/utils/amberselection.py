# strings = ['/A-E:34-45@CA,BB', '  /A-E :34-37,47-49,TRP,ACE @CA,CA,CA ', ':34-47']
# for string in strings: print(amberSelection(string))

import numpy as np
from   .protein_sel import three2one

def strsplit2dict(string, sep=['/', ':', '@']):
    string   = string.strip()
    data     = {}
    save     = ''
    prev_sep = None

    for s in string:
        if s in sep:
            if save.strip():
                data[prev_sep] = save.strip()
            save = ''
            prev_sep = s
            continue
        else:
            save += s
    data[prev_sep] = save.strip()
    return data

def dictsplit(dic):
    data = {}
    for key, value in dic.items():
        data[key] = []
        tmp1 = [s.strip() for s in value.split(',') if s.strip()]
        tmp2 = [s.split('-') for s in tmp1]
        for t in tmp2:
            if len(t) == 1:
                data[key].append(int(t[0]) if t[0].isdigit() else t[0])
            elif len(t) > 1.5:
                if t[0].isdigit():
                    # resid
                    tlist = np.arange(int(t[0]), int(t[-1]) + 1, dtype=np.int64)
                else:
                    # chain
                    tlist = np.arange(ord(t[0]), ord(t[-1]) + 1).astype('uint32').view('U1')
                data[key].extend(list(tlist))
    return data

def changekey(dic):
    data = {}
    for key, value in dic.items():
        if key == '/':
            data['chain'] = value
        elif key == ':':
            for val in value:
                if isinstance(val, str):
                    if 'resname' not in data.keys():
                        data['resname'] = [val]
                    else:
                        data['resname'].append(val)
                else:
                    if 'resid' not in data.keys():
                        data['resid'] = [int(val)]
                    else:
                        data['resid'].append(int(val))
        elif key == '@':
            data['name']  = value
    return data


def amberSelection(string):
    if string == 'protein':
        data = {'resname': list(three2one.keys())}
    else:
        data = changekey(dictsplit(strsplit2dict(string)))
    return data
    
def separateStarNoStar(value):
    valstar   = []
    valnostar = []
    for val in value:
        sval = val.split('*')
        if len(sval) > 1.5:
            valstar.append(sval[0])
        else:
            valnostar.append(val)
    return valstar, valnostar

def strsplitbysepinc(string, sep=['(', ')', '&', '|']):
    string = string.strip()
    data = []
    save = ''
    for s in string:
        if s in sep:
            if save.strip(): data.append(save.strip())
            if s.strip(): data.append(s.strip())
            save = ''
        else:
            save += s

    if save.strip(): data.append(save.strip())
    return data

def strsplitbysepexc(string, sep=['(', ')', '&', '|']):
    string = string.strip()
    data = []
    save = ''
    for s in string:
        if s in sep:
            if save.strip(): data.append(save.strip())
            save = ''
        else:
            save += s

    if save.strip(): data.append(save.strip())
    return data


