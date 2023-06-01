def pandas2dict(df, data):
	datakeys = data.keys()
	dfkeys   = df.columns
	ckeys    = set(datakeys).intersection(dfkeys)

	for ckey in ckeys:
		data[ckey].extend(df[ckey].tolist())

	return data
