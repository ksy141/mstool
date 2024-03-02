def saveargs(ifile, args, exclude=[]):
    if not isinstance(exclude, list):
        exclude = [exclude]

    with open(ifile, 'w') as W:
        for key, value in vars(args).items():
            if key not in exclude:
                W.write(f'{key:20} = {value}\n')
