def amberSelection(s):
    '''amber selection to dict.
    
    Parameters
    ----------
    s : str

    Returns
    -------
    d : dict

    Examples
    --------
    >>> amberSelection('/A:362@CA')
    >>> {'name': 'CA', 'chain': 'A', 'resid': 362}
    >>> amberSelection('/A:362')
    >>> {'name': None, 'chain': 'A', 'resid': 362}
    '''

    chain = None
    resid = None
    name  = None

    separators = ['/', ':', '@']
    collect = []
    tmp = ''
    for char in s:
        if char in separators:
            collect.append(tmp)
            tmp = ''

        tmp += char
    collect.append(tmp)
    collect.pop(0)

    for c in collect:
        if c.startswith('/'):
            chain = c[1:]
        elif c.startswith(':'):
            resid = int(c[1:])
        elif c.startswith('@'):
            name = c[1:]
        else:
            assert 0 == 1, 'unrecogonized separator, %s' %c[0]

    return {'name': name, 'chain': chain, 'resid': resid}


